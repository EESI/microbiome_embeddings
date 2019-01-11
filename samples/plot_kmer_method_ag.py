from sklearn.externals import joblib
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
k = 6
[X, y] = joblib.load('ag' + str(k) + '_mer_reduced_tsn.pkl')
# [X, y] = joblib.load('ag' + str(k) + '_mer_reduced_sub.pkl')
X /= X.std(axis=0)

font = {'size'   : 10}

matplotlib.rc('font', **font)
color_list = ['r', 'g', 'b']
legend_list = ['UBERON:feces', 'UBERON:tongue', 'UBERON:skin']
# X /= X.std(axis=0)
plt.figure(figsize = (7, 7))
for i in [0, 2, 1]:
    plt.scatter(X[y==i, 0], X[y==i, 1], color = color_list[i], label = legend_list[i], s = 5, alpha = 0.5)

lgnd = plt.legend(bbox_to_anchor=(0.5,-0.15), loc="lower center", ncol=3)
lgnd.legendHandles[0]._sizes = [60]
lgnd.legendHandles[1]._sizes = [60]
lgnd.legendHandles[2]._sizes = [60]
lgnd.legendHandles[0].set_alpha(1)
lgnd.legendHandles[1].set_alpha(1)
lgnd.legendHandles[2].set_alpha(1)
plt.grid()
# plt.title('ICA on ' + str(k) + '-mer frequency table')
plt.xlabel('Axis 1')
plt.ylabel('Axis 2')
plt.savefig('ag_tsn_' + str(k) + 'mr.pdf', format='pdf', dpi=1000, bbox_inches='tight')
plt.show()
