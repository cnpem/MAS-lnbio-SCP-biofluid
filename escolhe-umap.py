# Calculate sequence difference view for UMAP projections with given
# metrics, initializations and seeds.

# Usage example is:
# python escolhe-umap.py saliva.tsv euclidean,cosine random,pca 0 1000 | tee escolhe-umap.txt}

# Copyright 2025 Guilherme P. Telles.
# This program is distributed under the terms of WTFPL v2.

import sys

import pandas as pd
import numpy as np

from umap import UMAP

import warnings
warnings.filterwarnings("ignore")

from sklearn.neighbors import NearestNeighbors



tsvf = sys.argv[1]
metrics = sys.argv[2].split(',')
inits = sys.argv[3].split(',')
seedi = int(sys.argv[4])
seedf = int(sys.argv[5])


#fluid = tsvf.split('-')[0]
#print(tsvf,metrics,inits,seedi,seedf,fluid)


raw = pd.read_csv(tsvf, sep="\t")
raw = raw.drop(["PG.Genes","PG.FastaFiles"], axis=1)
raw = raw.fillna(0)

t = raw.set_index("PG.ProteinGroups").transpose()

D = t.values


'''
metrics = ['euclidean']
inits = ['random']
seedi = 1
seedf = 2

aligned = {'Dim1': [1, 2, 3, 4, 5],
           'Dim2': [6, 7, 8, 9, 10],
           'Dim3': [11, 12, 13, 14, 15]}
D = pd.DataFrame(aligned)
print(D)
'''


n,m = D.shape

nn = NearestNeighbors(metric='euclidean').fit(D)
distances, indices = nn.kneighbors(D)
originalNN = indices


k = n//2  # the k in sequence difference view equation


roN = list()
for i in range(n) :
    roN.append(list(originalNN[i][1:]))


for init in inits :

    for metric in metrics :

      for seed in range(seedi,seedf) :

        umap = UMAP(metric=metric, init=init, random_state=seed)
        umap_result = umap.fit_transform(D)
        
        P = np.column_stack((umap_result[:, 0], umap_result[:, 1]))

        nn = NearestNeighbors(metric='euclidean').fit(P)
        distances, indices = nn.kneighbors(P)
        projectedNN = indices

        ro2 = list()
        for i in range(n) :
            ro2.append(list(projectedNN[i][1:]))
            
        #print('ro2',ro2,sep='\n')
        #print('roN',roN,sep='\n')
        
        sumseqdv = 0.0
        
        for i in range(n) :

            s = 0.0
            for j in ro2[i][0:k] :
                s = s + (k - ro2[i].index(j)) * abs(ro2[i].index(j) - roN[i].index(j)) 

            sumseqdv += 0.5 * s
                
            s = 0.0
            for j in roN[i][0:k] :
                s = s + (k - roN[i].index(j)) * abs(ro2[i].index(j) - roN[i].index(j)) 

            sumseqdv += 0.5 * s

        '''
        sumsetdv = 0.0
        
        for i in range(n) :
            AN = list(originalNN[i][1:])
            A2 = list(projectedNN[i][1:])
        
            AN = AN[0:k]
            A2 = A2[0:k]
            
            I = list(set(A2) & set(AN))
            U = list(set(A2) | set(AN))
            
            sumsetdv += 1.0 / (len(I) / len(U))
        '''
        
        print(tsvf, metric, init, 'seed', seed, 'k', k, 'sequence_difference_view', sumseqdv)

exit()

