import glmnet_python
from glmnet import glmnet
import scipy, importlib, pprint, matplotlib.pyplot as plt, warnings
from glmnet import glmnet; from glmnetPlot import glmnetPlot
from glmnetPrint import glmnetPrint; from glmnetCoef import glmnetCoef; from glmnetPredict import glmnetPredict
from cvglmnet import cvglmnet; from cvglmnetCoef import cvglmnetCoef
from cvglmnetPlot import cvglmnetPlot; from cvglmnetPredict import cvglmnetPredict
import numpy as np
from sklearn.externals import joblib

k = 4

[train_kmer_table, train_kmer_label] = joblib.load('ag' + str(k)  + '_mer_train.pkl')
[test_kmer_table, test_kmer_label] = joblib.load('ag' + str(k)  + '_mer_test.pkl')
x = train_kmer_table/np.sum(train_kmer_table, axis = 1).reshape(-1,1)
y = train_kmer_label.copy()
#y = np.zeros((x.shape[0], 3))
#for i in range(x.shape[0]):
#    y[i, int(train_kmer_label[i])] = 1
print(x.shape, flush = True)
print(y.shape, flush = True)
del train_kmer_table
del train_kmer_label
x_test = test_kmer_table/np.sum(test_kmer_table, axis = 1).reshape(-1,1)
y_test = test_kmer_label.copy()
#y_test = np.zeros((x_test.shape[0], 3))
#for i in range(x_test.shape[0]):
#    y_test[i, int(test_kmer_label[i])] = 1
print(x_test.shape, flush = True)
print(y_test.shape, flush = True)
del test_kmer_table
del test_kmer_label

warnings.filterwarnings('ignore')
cvfit = cvglmnet(x = x.copy(), y = y.copy(), family = 'multinomial', ptype = 'class')
warnings.filterwarnings('default')
y_pred = cvglmnetPredict(cvfit, newx = x_test, s = 'lambda_min', ptype = 'class')

print("best alpha: ", cvfit['lambda_min'])
print(y_pred.shape, flush = True)
from sklearn.metrics import accuracy_score
print(accuracy_score(y_test, y_pred), flush = True)
