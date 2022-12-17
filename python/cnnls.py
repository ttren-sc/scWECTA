import torch
import random
import numpy as np

from sklearn import metrics

from scipy.optimize import nnls
from scipy.optimize import Bounds
from scipy.optimize import minimize
from scipy.optimize import fmin_slsqp
from scipy.optimize import least_squares

##  Non-negative least squares with equal sum to 1 constraint
class cnnls():
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.n = x.shape[1]

    ## Opitmization function
    def fun(self, w, x, y):
        pre = np.dot(x, w)
        loss = np.sum((pre - y) **2) / 2

        return loss

    ## Calculate weight matrix
    def weight(self):
        cons = list()
        w = np.random.random((self.n, 1))
        eq_cons = {
            "type": "eq",
            "fun": lambda w: np.sum(w) - 1
        }
        cons.append(eq_cons)

        bounds = Bounds(0, 1)
        res = minimize(self.fun, w, method="SLSQP", constraints=cons, options={'ftol': 1e-9, 'disp': False}, args=(self.x, self.y),
                       bounds=bounds)
        return res.x
