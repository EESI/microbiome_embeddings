import numpy as np
from sklearn.decomposition import fastica
from sklearn.externals import joblib

k = 10

[train_kmer_table, train_kmer_label] = joblib.load('ag' + str(k)  + '_mer_train.pkl')
[test_kmer_table, test_kmer_label] = joblib.load('ag' + str(k)  + '_mer_test.pkl')

X = np.concatenate((train_kmer_table, test_kmer_table), axis=0)
y = np.concatenate((train_kmer_label, test_kmer_label), axis=0)

X = X/np.sum(X, axis = 1).reshape(-1,1)
del train_kmer_table
del test_kmer_table

_, _, reduced = fastica(X, n_components=2, algorithm="deflation", whiten=True,
            fun="logcosh", fun_args=None, max_iter=200, tol=1e-04, w_init=None,
            random_state=3, return_X_mean=False, compute_sources=False,
            return_n_iter=False)

# can be removed
#reduced /= reduced.std(axis=0)

joblib.dump([reduced, y], 'ag' + str(k)  + '_mer_reduced.pkl')
