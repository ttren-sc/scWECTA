import pandas as pd
import random
from sklearn.model_selection import train_test_split


class findMarkers():
    def __init__(self, path, speciesType = "Human", tissueType = "Pancreas"):
        self.speciesType = speciesType
        self.tissueType = tissueType
        # self.markers = pd.read_table(path)
        self.path = path

    def _findMarkers(self):
        self._check_params()
        markers = pd.read_table(self.path)
        # Feature selection
        ## Get markers of pancreas
        tissue_marker = markers.loc[(markers["speciesType"] == self.speciesType) & (markers["tissueType"] == self.tissueType), 'geneSymbol']
        # tissue_marker = self.markers.loc[self.markers[self.tissueType] == self.tissueType, 'geneSymbol']
        tissue_marker.index = range(len(tissue_marker))
        marker_genes = []
        for i in range(len(tissue_marker)):
            if pd.isna(tissue_marker[i]):
               continue
            marker = str.split(tissue_marker[i], ", ")
            marker_genes.append(marker)
        marker_genes = list(set([str(x) for item in marker_genes for x in item])) # 115, 101

        return marker_genes

    def _check_params(self):
        markers = pd.read_table(self.path)
        speciesList = list(set(markers.loc[:, 'speciesType']))  # 2: Human, Mouse
        tissueList = list(set(markers.loc[:, 'tissueType']))  # 119


        if self.speciesType not in speciesList:
            raise Exception('Species is not available, please input species from the following list: {}.'.format(speciesList))
        if self.tissueType not in tissueList:
            raise Exception('Tissue is not available, please input tissue from the following list: {}.'.format(tissueList))

# gene selection by padding 0 to test
class geneSelection():
    def __init__(self, train, validation, test, featureSet):
        # self.type = type
        self.train = train
        self.validation = validation
        self.test = test
        self.select_features = featureSet

    def feature_filter(self):
        train_features = list(set(self.train.columns).intersection(self.select_features))
        return train_features

    def _reSort(self, df):
        return df.reindex(sorted(df.columns), axis=1)


    def filtered_result(self):
        # self.type is True when featureSet is marker genes
        # if self.type:
        #     train_features = self.feature_filter()
        #     print(train_features)
        # else:
        #     train_features = self.select_features
        train_features = self.feature_filter()
        test_features = self.test.columns
        intersect = list(set(test_features).intersection(train_features))
        miss = list(set(train_features).difference(intersect))

        if not miss:
            return self.train.loc[:, train_features], self.validation.loc[:, train_features], self.test.loc[:, train_features] # self.train.loc[:, train_features],
        # else:
        #     result = self.test.loc[:, intersect]
        #     for item in miss:
        #         result[miss] = 0
        #     return self._reSort(self.train.loc[:, train_features]), self._reSort(result) # self._reSort(self.train.loc[:, train_features]),
        else:
            # result = self.test.loc[:, intersect]
            test = pd.DataFrame(index=self.test.index, columns=train_features)
            for i in train_features:
                if i in miss:
                    test[i] = 0
                else:
                    test[i] = self.test[i]
            return self.train.loc[:, train_features], self.validation.loc[:, train_features], test

# gene selection by intersection
class geneSelection2():
    def __init__(self, train, validation, test, featureSet):
        self.train = train
        self.validation = validation
        self.test = test
        self.select_features = featureSet

    def _reSort(self, df):
        return df.reindex(sorted(df.columns), axis=1)

    def filtered_result(self):
        features = list(set(self.train.columns).intersection(self.select_features))

        return self.train.loc[:, features], self.validation.loc[:, features], self.test.loc[:, features]
