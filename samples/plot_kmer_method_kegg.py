f = open('seqs_tax_assignments.txt', 'rb')
tax2label = {}
tax_dict = {}
tmp = 0
unc = 0
for line in f:
    med = line.decode('utf-8').split('\t')
    tax_labels = med[1].split(';')
    for item in tax_labels:
        level_tax = item.split('__')
        if len(level_tax) == 2:
            level = level_tax[0]
            tax = level_tax[1]
            if level in tax2label:
                if tax not in tax2label[level]:
                    tax2label[level][tax] = tmp
                    tmp += 1
            else:
                tax2label[level] = {}
                tax2label[level][tax] = tmp
                tmp += 1
        else:
            unc += 1
    tax_dict[med[0]] = tax_labels
f.close()

import numpy as np
from sklearn.externals import joblib
import matplotlib.pyplot as plt
import matplotlib
k = 6
X = np.load('kegg_' + str(k) + '_mer_reduced_tsn.npy')
# X = joblib.load('kegg_' + str(k) + '_mer_reduced.pkl')
y = joblib.load('kegg_10_mer_ids.pkl')

intest = ['Actinobacteria', 'Bacteroidetes', 'Chlamydiae', 'Cyanobacteria', 'Euryarchaeota', 'Firmicutes', 'Proteobacteria','Tenericutes']
y_text = []
for item in y:
    try:
        tmp = tax_dict[item][1]
        tag, name = tmp.split('__')
        if tag == 'p':
            if name in intest:
                y_text.append(tax2label['p'][name])
            else:
                y_text.append(-1)
        else:
            y_text.append(-1)
    except:
        y_text.append(-1)
        
c_idx = []
for item in intest:
    c_idx.append(tax2label['p'][item])

y = np.array(y_text)
X = X[y!=-1, :]
y = y[y!=-1]


font = {'size'   : 10}

matplotlib.rc('font', **font)
# color_list = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'w']
legend_list = intest
# X /= X.std(axis=0)
plt.figure(figsize = (7, 7))
for idx, i in enumerate(c_idx):
    plt.scatter(X[y==i, 0], X[y==i, 1], label = legend_list[idx], s = 5, alpha = 0.5)

lgnd = plt.legend(bbox_to_anchor=(1.15,0.5), loc="center")
for item in lgnd.legendHandles:
    item._sizes = [60]
    item.set_alpha(1)

plt.grid()
# plt.title('ICA on ' + str(k) + '-mer frequency table')
plt.xlabel('Axis 1')
plt.ylabel('Axis 2')
plt.savefig('kegg_tsn_' + str(k) + 'mr.pdf', format='pdf', dpi=1000, bbox_inches='tight')
plt.show()


interst_family = ['Bacillaceae', 'Enterobacteriaceae', 'Streptococcaceae']
for ix, interst_name in enumerate(interst_family):
    X = np.load('kegg_' + str(k) + '_mer_reduced_tsn.npy')
    y = joblib.load('kegg_10_mer_ids.pkl')
    y_text = []
    for item in y:
        try:
            tmp = tax_dict[item][4]
            tag, name = tmp.split('__')
            if tag == 'f':
                if name == interst_name:
                    tmp = tax_dict[item][5]
                    tag, name = tmp.split('__')
                    if tag == 'g':
                        y_text.append(tax2label['g'][name])
                    else:
                        y_text.append(-1)
                else:
                    y_text.append(-1)
            else:
                y_text.append(-1)
        except:
            y_text.append(-1)
    y = np.array(y_text)
    X = X[y!=-1, :]
    y = y[y!=-1]
    c_idx = np.unique(y)
    
    font = {'size'   : 20}

    matplotlib.rc('font', **font)

    plt.figure(figsize = (10, 5))
    for idx, i in enumerate(c_idx):
        plt.scatter(X[y==i, 0], X[y==i, 1], s = 15, alpha = 0.5)

    plt.grid()
    # plt.title('ICA on ' + str(k) + '-mer frequency table')
    plt.xlabel('Axis 1')
    plt.ylabel('Axis 2')
    if ix == 0:
        plt.xlim(-60, 10)
    if ix == 1:
        plt.xlim(-35, 35)
    if ix == 2:
        plt.xlim(-20, 50)
    plt.savefig('kegg_tsn_' + interst_name + '_' + str(k) + 'mr.pdf', format='pdf', dpi=1000, bbox_inches='tight')
    plt.show()
