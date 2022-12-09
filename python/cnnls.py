import torch
import random
import numpy as np

from sklearn import metrics

from scipy.optimize import nnls
from scipy.optimize import Bounds
from scipy.optimize import minimize
from scipy.optimize import fmin_slsqp
from scipy.optimize import least_squares


##  Set random seed
def init_random_seed(manual_seed):
    seed = None
    if manual_seed is None:
        seed = random.randint(1,10000)
    else:
        seed = manual_seed
    print("use random seed: {}".format(seed))
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

##  Non-negative least squares with equal sum to 1 constraint
class cnnls():
    def __init__(self, x, y, seed = 1):
        self.x = x
        self.y = y
        self.n = x.shape[1]
        self.seed = seed

    ## Set seed
    def init_random_seed(self):
        seed = None
        if self.seed is None:
            seed = random.randint(1,10000)
        else:
            seed = self.seed
        print("use random seed: {}".format(seed))
        random.seed(seed)
        np.random.seed(seed)
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)

    ## Opitmization function
    def fun(self, w, x, y):
        pre = np.dot(x, w)
        loss = np.sum((pre - y) **2) / 2

        return loss

    ## Calculate weight matrix
    def weight(self):
        self.init_random_seed()
        cons = list()
        w = np.random.random((self.n, 1))
        eq_cons = {
            "type": "eq",
            "fun": lambda w: np.sum(w) - 1
        }
        cons.append(eq_cons)

        # e = 1e-10
        # ineq_cons = {
        #     "type": "ineq",
        #     "fun": lambda w: w
        # }
        # cons.append(ineq_cons)

        bounds = Bounds(0, 1)
        res = minimize(self.fun, w, method="SLSQP", constraints=cons, options={'ftol': 1e-9, 'disp': True}, args=(self.x, self.y),
                       bounds=bounds)
        return res.x


## Test with random data
if __name__ == '__main__':
    init_random_seed(1)
    x = np.random.random((10, 5))
    y = np.random.randint(0, 2, 10)

    w = cnnls(x, y, 1).weight()
    print(w)
