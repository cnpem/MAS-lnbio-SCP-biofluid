# Print the UMAP projection of saliva or tear with given metric,
# initialization and seed.

# Also creata a plot for each protein in saliva-fda.tsv or
# tear-fda.tsv, depcting intensity as color saturation.

# Usage example is:
# python plot.py saliva.tsv euclidean random 880
# which will read saliva-fda.tsv also.

# Copyright 2025 Guilherme P. Telles.
# This program is distributed under the terms of WTFPL v2.


import re
import sys

import numpy as np

import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1 import make_axes_locatable

from umap import UMAP

import warnings
warnings.filterwarnings("ignore")


tsvf = sys.argv[1]
metric = sys.argv[2]
init = sys.argv[3]
seed = int(sys.argv[4])

prefix = tsvf.replace('.tsv','')
fluid = prefix.split('-',1)[0]

#print(tsvf, metric, init, seed, prefix, fluid)



# Load the size of saliva cells:
df = pd.read_csv("saliva-size.tsv", sep="\t")

ssizes = {}
for index,row in df.iterrows() :
  ssizes[row['Position_Plate']] = row['Diameter (uM)']

#print(ssizes)


# Load the size of tear cells:
df = pd.read_csv("tear-size.tsv", sep="\t")

tsizes = {}
for index,row in df.iterrows() :
  tsizes[row['Position_Plate']] = row['Diameter (uM)']

#print(tsizes)
condition_colors = {
    "Tears" : "#3466AC",
    "Saliva" : "#D9762B"
}



data = pd.read_csv(tsvf, sep="\t")
data = data.drop(["PG.Genes","PG.FastaFiles"], axis=1)
data = data.fillna(0)
data = data.set_index("PG.ProteinGroups").transpose()



cells = []
conditions = []
colors = []
sizes = []

for r in data.index :

    cond = re.findall("Biofluid_([a-zA-Z]+)_", r)[0]
    conditions.append(cond)

    cell = re.findall("_([a-zA-Z0-9]+)\.", r)[0]
    cells.append(cell)

    colors.append(condition_colors[cond])

    if cond == 'Saliva' :
        sizes.append(3.14*ssizes[cell]*ssizes[cell]/4)
    else :
        sizes.append(3.1416*tsizes[cell]*tsizes[cell]/4)

# print(cells)
# print(conditions)
# print(colors)
# print(sizes)

# umap:
umap = UMAP(metric=metric, init=init, random_state=seed)
umap_result = umap.fit_transform(data.values)

plt.figure(figsize=(10,8))
sc = plt.scatter(umap_result[:, 0], umap_result[:, 1], edgecolors='black', s=sizes, c=colors)

#bar = plt.colorbar(sc,fraction=0.05,shrink=.5,pad=0.05,aspect=15)
#bar.set_label("Intensity")  

#for i, txt in enumerate(cells):
#  plt.annotate(txt, (umap_result[:,0][i],umap_result[:,1][i]))
  
plt.title("UMAP Projection of " + prefix)
plt.xlabel("UMAP Component 1")
plt.ylabel("UMAP Component 2")

  
'''
legend_handles = [
mpatches.Patch(color=color, label=condition)
for condition, color in condition_colors.items()
]

plt.legend(handles=legend_handles, title="Conditions")
'''
  
#plt.show()

plt.savefig(prefix + '-' + metric + '-' + init + '-' + str(seed) + \
            '.svg',format='svg') #bbox_inches='tight'


biof = fluid + "-fda.tsv"

bio = pd.read_csv(biof,sep="\t")
bio = bio.fillna(0)


n = len(bio)

for i in range(0,n) :
  z = list(bio.iloc[i, 0:])
  prot = z[0]
  ints = z[3:]

  print(prot)

  m = min(ints)
  M = max(ints)
  
  for j in range(0,len(ints)) :
    ints[j] = (ints[j]-m)/M

  if cond == 'Saliva' :
    cm = plt.colormaps.get_cmap('Oranges')
  else :
    cm = plt.colormaps.get_cmap('Blues')
  
  plt.figure(figsize=(10, 8))
  sc = plt.scatter(umap_result[:, 0], umap_result[:, 1], edgecolors='black', c=ints, s=sizes, cmap=cm)

  bar = plt.colorbar(sc,fraction=0.05,shrink=.5,pad=0.05,aspect=15)
  bar.set_label("Intensity")  

  # for i, txt in enumerate(cells):
  #   plt.annotate(txt, (umap_result[:,0][i],umap_result[:,1][i]))
    
  plt.title("UMAP Projection of " + prefix + " colored by " + prot)
  plt.xlabel("UMAP Component 1")
  plt.ylabel("UMAP Component 2")

  
  '''
  legend_handles = [
    mpatches.Patch(color=color, label=condition)
      for condition, color in condition_colors.items()
  ]
  
  plt.legend(handles=legend_handles, title="Conditions")
  '''
  
  # plt.show()

  plt.savefig(prefix + '-' + prot + '-' + metric + '-' + init + '-' +\
              str(seed) + '.svg',format='svg') #bbox_inches='tight'
