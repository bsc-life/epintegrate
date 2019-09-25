import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
from sklearn.cluster import KMeans
from scipy.stats.stats import pearsonr
import itertools

df = pd.read_csv('complete_genome_data.csv', sep='\t')

cells = ['CLL_12', 'CLL_110', 'CLL_1525', 'NBC', 'MBC', 'GCBC', 'PBC']
num_classes = len(cells)
genome_length = len(df[(df.Cell == 'CLL_12') & (df.Repl == 12)]) #canviar-ho per fer-ho general

def calinski_harabasz_score(X, labels):
        n_samples, _ = X.shape
        n_labels = len(set(labels))

        extra_disp, intra_disp = 0., 0.
        mean = np.mean(X, axis=0)
        for k in range(n_labels):
                cluster_k = X[labels == k]
                mean_k = np.mean(cluster_k, axis=0)
                extra_disp += len(cluster_k) * np.sum((mean_k - mean) ** 2)
                intra_disp += np.sum((cluster_k - mean_k) ** 2)

        return (1. if intra_disp == 0. else extra_disp * (n_samples - n_labels) / (intra_disp * (n_labels - 1.)))

def get_importance_vector(cl1,cl2):
    importance_vector = np.zeros(genome_length)
    occurrence_vector = np.zeros(genome_length)
    num_samples = 100000 #100000
    remaining_percentage = 0.01 #We only take 1% of the data. And the clusters still form!
    remaining_number = int(genome_length*remaining_percentage)

    
    for i in range(num_samples):
        #len-1 because length not coinciding in all replicates
        remaining_positions = random.sample(range(genome_length), remaining_number)
        xx = []

        for i in range(len(remaining_positions)):
            occurrence_vector[remaining_positions[i]] += 1
            
        for c in [cl1, cl2]:
            for r in range(1,21): #canviar per unique values in column Repl per fer-ho general
                #print(c)
                #print(r)
                dftmp = df[(df.Cell == c) & (df.Repl == r)]
                dftmp = dftmp.reset_index()
                #print(len(dftmp))
                #v = [dftmp.loc[j]['Value'] for j in remaining_positions] #què passa si no existeix? com afegir Nan?
                v = []
                for j in remaining_positions:
                    if j < len(dftmp):
                        #print(dftmp.loc[j]['Value'])
                        v.append(dftmp.loc[j]['Value'])
                    else:
                        #print("No")
                        v.append(0.5)
                xx.append(v)

        xx = np.asarray(xx)

        kmeans = KMeans(n_clusters=2)
        kmeans.fit(xx)
        
        cluster_performance = calinski_harabasz_score(xx, kmeans.labels_)

        for index in remaining_positions:
            importance_vector[index] += cluster_performance

    for i in range(len(importance_vector)):
        if (occurrence_vector[i] != 0):
            importance_vector[i] = importance_vector[i] / occurrence_vector[i]
        else:
            importance_vector[i] = 0
    #print(importance_vector)
    return importance_vector

imp_matrix = np.zeros((num_classes, num_classes, genome_length))

for i in range(num_classes):
        for j in range(i+1, num_classes):
                imp_matrix[i][j] = get_importance_vector(cells[i], cells[j])
                filename = 'iv_' + cells[i] + '_' + cells[j] + '.txt'
                np.savetxt(filename, imp_matrix[i][j])

#np.savetxt(filename, imp_matrix)
