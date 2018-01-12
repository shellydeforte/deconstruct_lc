import matplotlib.pyplot as plt
import numpy as np
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
    score = clf.score(X, y)
    print(score)
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
    score = clf.score(X, y)
    print(score)
    return clf


def plot_contour(clf, X, y):
    plt.scatter(X[:, 0], X[:, 1], c=y, zorder=9, cmap=plt.cm.Paired,
                edgecolor='black', s=5)
    x_min = X[:, 0].min()  # x and y are each of the classes
    x_max = X[:, 0].max()
    y_min = X[:, 1].min()
    y_max = X[:, 1].max()
    # Creates a grid with the same number on each line
    XX, YY = np.mgrid[x_min:x_max:200j, y_min:y_max:200j]
    # XX, YY are the individual points of the grid.
    # np.c_[XX.ravel(), YY.ravel()] creates the x,y pairs of the grid (I think)
    Z = clf.decision_function(np.c_[XX.ravel(), YY.ravel()])
    # Put the result into a color plot
    Z = Z.reshape(XX.shape)
    plt.pcolormesh(XX, YY, Z > 0, cmap=plt.cm.Paired)
    plt.contour(XX, YY, Z, colors=['white'],
                linestyles=['-'], linewidth=5, levels=[0], zorder=10)
    #plt.xlim([100, 500])
    #plt.ylim([-20, 100])
    plt.xlabel('LCA score')
    plt.ylabel('LCE score')
    plt.show()


def main():
    pass


if __name__ == '__main__':
    main()
