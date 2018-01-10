import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.svm import SVC
from sklearn.svm import LinearSVC

def smooth_rbf(X, y):
    clf = SVC(C=0.1, cache_size=200, class_weight=None, coef0=0.0,
              decision_function_shape='ovr', degree=3, gamma='auto',
              kernel='rbf', max_iter=-1, probability=False, random_state=None,
              shrinking=True, tol=0.001, verbose=False)
    clf.fit(X, y)
    score = clf.score(X, y)
    return clf

def linear_svc(X, y):
    clf = LinearSVC(random_state=0)
    clf.fit(X, y)
    return clf.score(X, y)