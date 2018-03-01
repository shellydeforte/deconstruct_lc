"""
Created by Shelly DeForte, Michnick Lab, University of Montreal 2017-2018
https://github.com/shellydeforte/deconstruct_lc/
"""
from sklearn.svm import SVC


def smooth_rbf(X, y):
    clf = SVC(kernel='rbf',
              C=0.1,
              cache_size=500,
              class_weight=None,
              random_state=0,
              decision_function_shape='ovr',
              gamma='auto',
              max_iter=-1,
              probability=False,
              shrinking=True,
              tol=0.001,
              verbose=False)
    clf.fit(X, y)
    return clf


def linear_svc(X, y):
    clf = SVC(kernel='linear',
              C=1,
              cache_size=500,
              class_weight=None,
              random_state=None,
              decision_function_shape='ovr',
              gamma='auto',
              max_iter=-1,
              probability=False,
              shrinking=True,
              tol=.001,
              verbose=False)
    clf.fit(X, y)
    return clf