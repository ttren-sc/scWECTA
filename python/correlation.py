import time
import math
import random
import numpy as np
import pandas as pd
import scipy.stats
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import StackingClassifier
from sklearn.ensemble import VotingClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.metrics import accuracy_score
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted
from sklearn.utils import column_or_1d
from sklearn.preprocessing import LabelEncoder
from sklearn.utils import indexable
import warnings

warnings.filterwarnings('ignore')

class Correlation(BaseEstimator, ClassifierMixin):
    def __init__(self, similarity='spearman'):
        self.similarity = similarity

    def fit(self, X, y=None):
        self._check_params(X, y)
        X, y = check_X_y(X, y)
        X, y = indexable(X, y)

        self.X = X

        classes = list(set(y))
        classes.sort()
        self.classes_ = classes
        self.encode = LabelEncoder()
        self.encode.fit(classes)
        class_encode = self.encode.transform(y)
        self.y = class_encode

        self.df = pd.DataFrame(y)
        self.df.columns = ['class']

        return self

    def predict(self, X):
        check_is_fitted(self)
        X = check_array(X)
        y_pred = np.argmax(self.predict_proba(X), axis=1)
        result = self.encode.inverse_transform(y_pred)

        return result

    def predict_proba(self, X):
        check_is_fitted(self)
        X = check_array(X)
        self.pred = X

        ref = pd.DataFrame(self.X)
        reference = pd.DataFrame(index=range(ref.shape[1]), columns=self.classes_)

        for i in self.classes_:
            cells = self.df.loc[self.df['class'] == i].index.tolist()
            reference.loc[:, i] = np.array(ref.loc[cells, :].sum(axis=0)/ len(cells))
        reference.columns = range(len(self.classes_))

        y_prob = None

        if self.similarity == 'cosine':
            y_prob = cosine_similarity(np.array(reference).T, X)
            y_prob = y_prob.T

        elif self.similarity == "pearson":
            y_prob = reference.apply(self.pearson, 0)

        elif self.similarity == 'spearman':
            y_prob = reference.apply(self.spearman, 0)

        y_prob = np.array(y_prob)

        return y_prob


    def pearson(self, x):
        result = np.apply_along_axis(scipy.stats.pearsonr, 1, self.pred, np.array(x))

        return result[:, 0]

    def spearman(self, x):
        result = np.apply_along_axis(scipy.stats.spearmanr, 1, self.pred, np.array(x))

        return result[:, 0]

    def _check_params(self, X, y):
        if self.similarity not in ['cosine', 'pearson', 'spearman']:
            raise Exception('"similarity" should be chosen among "Cosine", "Pearson", or "Spearman"')


