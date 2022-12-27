import os
import re
import time
import math
import random

import numpy as np
import pandas as pd

import scipy.stats
from sklearn.utils.validation import check_is_fitted

from sklearn.preprocessing import LabelEncoder

from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LinearRegression

from sklearn.ensemble import StackingClassifier
from sklearn.ensemble import VotingClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.decomposition import NMF
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.model_selection import KFold

from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score, auc, f1_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import cohen_kappa_score

from correlation import Correlation
from featureSelection import geneSelection
from featureSelection import findMarkers
from cnnls import cnnls

from collections import Counter
import argparse
import warnings

probs = ['lr', 'lda', 'nb', 'ann', 'svm'] # Models trained under a probabilistic framework.
warnings.filterwarnings('ignore')
#######------- Run all the steps once with train dataset splitted into two datasets and threshold -------#######

class Ensemble(BaseEstimator, ClassifierMixin):
    def __init__(
            self,
            estimators,
            geneSets=None,
            threshold = 0.5,
            calibrate = "sigmoid"
    ):
        # super().__init__(estimators=estimators)
        self.estimators = estimators
        self.geneSets = geneSets
        self.threshold = threshold
        self.calibrate = calibrate

    def fit(self, Xd, yd, sample_weight=None):
        # self._check_params(X, y)
        skf = StratifiedKFold(n_splits=2)
        skf.get_n_splits(Xd, yd)
        for train_index, test_index in skf.split(Xd, yd):
            self.X_train, self.X_test = Xd.iloc[train_index], Xd.iloc[test_index]
            self.y_train, self.y_test = yd.iloc[train_index], yd.iloc[test_index]

        self.df = yd

        # self.xdim = len(self.X_[0])
        self.n = len(yd['label'])

        classes = list(set(yd['label']))
        classes.sort()
        self.classes_ = classes
        self.encode = LabelEncoder()
        self.encode.fit(classes)
        # class_encode = self.encode.transform(yd)
        self.y_train_label = self.encode.transform(self.y_train['label'])
        self.y_test_label = self.encode.transform(self.y_test['label'])

        self.fitted_ = True

        return self

    def predict(self, X):
        check_is_fitted(self)
        probas, weight_ = self.predict_proba(X) # probas is the matrix stored the predicted probability of every cell type for every cell, weight_ is the cell-type specific weight matrix.
        pred_ensem = np.argmax(probas, axis=1)
        pred_ensem_max = np.max(probas, axis=1)

        # result = np.array([self.classes_[i] for i in pred_ensem])
        result = self.encode.inverse_transform(pred_ensem)

        for i in range(len(pred_ensem)):
            if pred_ensem_max[i] < self.threshold:
                result[i] = "unassigned"

        return result

    def _collect_probas(self, X):
        """Collect results from clf.predict calls."""
        meta_train = []
        meta_test = []
        count = 0
        for name, model in self.estimators.items():

            print("Now is classifier: {}".format(count))
            genes = self.geneSets[count]
            # train = self.X_.loc[:, genes]
            train_train, train_test, test = geneSelection(train=self.X_train, validation = self.X_test, test=X, featureSet=genes).filtered_result()

            # test = np.array(test)
            if model in probs:
                model.fit(train_train, self.y_train_label)
                y_prob1 = model.predict_proba(train_test)
                y_prob2 = model.predict_proba(test)
                meta_train.append(y_prob1)
                meta_test.append(y_prob2)
            else:
                # model.fit(self.X_, self.y_)
                model = CalibratedClassifierCV(model, method=self.calibrate, cv=2)  # sigmoid, isotonic
                model.fit(train_train, self.y_train_label)
                y_prob1 = model.predict_proba(train_test)
                y_prob2 = model.predict_proba(test)
                meta_train.append(y_prob1)
                meta_test.append(y_prob2)
            count += 1
        return meta_train, meta_test     # np.asarray([clf.predict_proba(X) for clf in self.estimators_])

    def predict_proba(self, X):
        check_is_fitted(self)

        x_train, x_test = self._collect_probas(X)
        n = len(self.estimators)
        gt = pd.DataFrame(0, index=self.X_test.index, columns=self.classes_)
        for i in self.classes_:
            cells = self.y_test.loc[self.y_test['label'] == i].index.tolist()
            gt.loc[cells, i] = 1
        k = len(self.classes_)
        weight = np.zeros((n, k))  # n, number of classifiers
        gt = np.array(gt)

        # b = np.zeros(k)
        ensemble = np.zeros((x_test[0].shape[0], k))
        # Predict cell-type specific weight
        for i in range(k):
            data = x_train[0][:, i]
            for j in range(1, n):
                data = np.vstack((data, x_train[j][:, i]))
            data = data.T
            # reg = LinearRegression(fit_intercept=True, positive=True).fit(data, gt[:, i])
            # w = reg.coef_
            # b[i] = reg.intercept_
            w = cnnls(data, gt[:, i]).weight()  # nnls with equal to 1 constraint
            weight[:, i] = w
        # Get soft ensemble of test dataset
        for i in range(k):
            tmp = 0
            for j in range(n):
                tmp += x_test[j][:, i] * weight[j, i]
            ensemble[:, i] += tmp
        weight_ = pd.DataFrame(weight)
        weight_.columns = self.classes_
        weight_.index = self.estimators
        return ensemble, weight_

## Get ensemble classifiers
def get_models(sim = "spearman"):
    models = dict()
    models['lr'] = LogisticRegression(solver="saga", penalty="elasticnet", l1_ratio=0.1, class_weight='balanced', max_iter=20)
    models['svm'] = SVC(class_weight='balanced', probability=True)
    models['knn'] = KNeighborsClassifier(weights="distance", n_neighbors=3)
    models['cart'] = DecisionTreeClassifier(class_weight='balanced')
    models['corr'] = Correlation(similarity=sim)
    return models

def get_ensemble(geneSets=None, threshold=0.5, sim = "spearman", calibrate = "sigmoid"):
    models = get_models(sim)
    model = Ensemble(estimators=models, geneSets= geneSets, threshold = threshold, calibrate = calibrate)
    return model

def create_dir_not_exist(path):
    if not os.path.exists(path):
        os.mkdir(path)

if __name__ == '__main__':
    # Input from linux.
    parser = argparse.ArgumentParser(description='scWECTA.')
    parser.add_argument('--train', help='The path of data relevant to Reference.')
    parser.add_argument('--test', help='The path of test data.')
    parser.add_argument('--marker', help='The path of marker gene list, and user need to make sure the file, called all_cell_markers.txt, is in this path')
    parser.add_argument('--output', help='The path used to store cell type annotation result.')
    parser.add_argument('-s', '--species', dest='species', default='Human', type=str, help='The name of speices, Human or Mouse')
    parser.add_argument('-t', '--tissue', dest='tissue', default='Pancreas', type=str, help='The name of tissue, tissue provided in CellMarker')
    parser.add_argument('--thred', default=0.5, type=float, help='Threshold for deciding final cell type annotation result, default is 0.5')
    parser.add_argument('--sim', default='spearman', type=str, help='Similarity metrics, including spearman, pearson and cosine.')
    parser.add_argument('--cal', default='sigmoid', type=str, help='Calibration method, including sigmoid and isotonic.')

    args = parser.parse_args()  # sys.argv[1:]

    train_path = args.train
    test_path = args.test
    marker_path = args.marker
    result_path = args.output
    speices = args.species
    tissue = args.tissue
    threshold = args.thred
    sim = args.sim
    calibrate = args.cal

    # Reference data.
    print("Reference data is loading.")

    train_dirs = os.listdir(train_path)
    train = list(filter(lambda x: re.match("train.csv", x) != None, train_dirs))
    Reference = np.loadtxt(fname=train_path + train[0], delimiter=',', dtype=pd.DataFrame)
    Reference = pd.DataFrame(Reference)
    genes = Reference.iloc[:, 0]
    genes = genes[1:]
    genes = genes.tolist()
    genes = [eval(i) for i in genes]  # remove the one of two quotation mark of each gene
    cells = Reference.iloc[0, :]
    cells = cells[1:]
    cells = cells.tolist()
    cells = [eval(i) for i in cells]
    Reference = Reference.iloc[1:, 1:]
    Reference.index = genes
    Reference.columns = cells

    # Label of training data.
    print("Label of reference data is loading.")

    ref_label = list(filter(lambda x: re.match("train_label.csv", x) != None, train_dirs))

    train_label = np.loadtxt(fname=train_path + ref_label[0], delimiter=',', dtype=pd.DataFrame)
    train_label = pd.DataFrame(train_label)
    cells = train_label.iloc[:, 0]
    cells = cells[1:]
    cells = cells.tolist()
    cells = [eval(i) for i in cells]  # remove the one of two quotation mark of each gene
    train_label = train_label.iloc[1:, 1:]
    ## make sure there is a column, named as 'label', in the label matrix of training data.
    train_label.columns = ['label']
    train_label.index = cells

    # Gene sets
    print("Gene sets of reference data are loading.")

    geneSets_file = list(set(train_dirs) - (set(train) | set(ref_label)))
    geneSets = list()
    for i in range(len(geneSets_file)):
        geneset = np.loadtxt(fname=train_path + geneSets_file[i], delimiter=',', dtype=pd.DataFrame)
        geneset = pd.DataFrame(geneset)
        geneSets.append(list(geneset.loc[:,1]))

    # Query data
    print("Query data is loading.")

    test_dirs = os.listdir(test_path)
    Query = np.loadtxt(fname=test_path + test_dirs[0], delimiter=',', dtype=pd.DataFrame)
    Query = pd.DataFrame(Query)
    genes = Query.iloc[:, 0]
    genes = genes[1:]
    genes = genes.tolist()
    genes = [eval(i) for i in genes]  # remove the one of two quotation mark of each gene
    cells = Query.iloc[0, :]
    cells = cells[1:]
    cells = cells.tolist()
    cells = [eval(i) for i in cells]
    Query = Query.iloc[1:, 1:]
    Query.index = genes
    Query.columns = cells

    print("Marker genes are loading.")

    marker_file = os.listdir(marker_path)
    markerGenes = findMarkers(speciesType=speices, tissueType=tissue, path= marker_path + marker_file[0])._findMarkers()
    geneSets.append(markerGenes)

    ## Random assign gene sets to different classifier
    random.shuffle(geneSets)
    count = []
    for i in range(5):
        count.append(len(geneSets[i]))

    ## Model training & test
    print("scWECTA is training.")

    model = get_ensemble(geneSets=geneSets, threshold=threshold, sim = sim, calibrate=calibrate)
    model.fit(Reference, train_label)
    # pred_prob, pred_ensem, pred_weight = model.predict(Query)
    pred_ensem = model.predict(Query)

    ## Predicting results
    create_dir_not_exist(result_path)
    df = pd.DataFrame({"label":pred_ensem})
    print(df)
    df.to_csv(result_path + "/anno.csv", index=False)
    print("scWECTA is over.")


