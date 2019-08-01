import os
import pandas as pd
import numpy as np
import random
from sklearn.cluster import KMeans
import matplotlib
import scipy.cluster.hierarchy as sch


data1 = []
data2 = []
data3 = []
data4 = []
data5 = []
data6 = []

experiments = ["S0159LH1", "S00W1BH1", "S00X9SH1", "S00XAQH1", "S00VKEH1", "S00Y8QH1", "S013ARH1", "S00W0DH1", "S00Y9OH1", "S015BHH1", "S015CFH1"]

labels = np.array([0,0,0,0,1,1,2,2,2,3,3])

def check_assignments(assignments):
    if (assignments[0] == assignments[1] == assignments[2] == assignments[3] and assignments[4] == assignments[5] and assignments[6] == assignments[7] == assignments[8] and assignments[9] == assignments[10]):
        return True
    else:
        return False

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

for exp in experiments:
	for file in os.listdir('histone_marks_data/'):
		if (file.startswith(exp) and 'chr1_' in file):
			fh = open('histone_marks_data/' + file)
			next(fh)
			next(fh)
			new1 = []
			new2 = []
			new3 = []
			new4 = []
			new5 = []
			new6 = []
			for line in fh:
				new1.append(int(line.split()[0]))
				new2.append(int(line.split()[1]))
				new3.append(int(line.split()[2]))
				new4.append(int(line.split()[3]))
				new5.append(int(line.split()[4]))
				new6.append(int(line.split()[5]))
			data1.append(new1)
			data2.append(new2)
			data3.append(new3)
			data4.append(new4)
			data5.append(new5)
			data6.append(new6)

for window_length in range(3, 200, 2):
	print(window_length)
	new_data1 = []
	new_data2 = []
	new_data3 = []
	new_data4 = []
	new_data5 = []
	new_data6 = []
    
	for k in range(len(data1)):
		new_data1.append([])
		new_data2.append([])
		new_data3.append([])
		new_data4.append([])
		new_data5.append([])
		new_data6.append([])

	for i in range(0, len(data1[0])-window_length-1, window_length):
		for j in range(len(data1)):
			pattern1 = data1[j][i:i+window_length]
			pattern2 = data2[j][i:i+window_length]
			pattern3 = data3[j][i:i+window_length]
			pattern4 = data4[j][i:i+window_length]
			pattern5 = data5[j][i:i+window_length]
			pattern6 = data6[j][i:i+window_length]

			suma1 = np.sum(pattern1)
			suma2 = np.sum(pattern2)
			suma3 = np.sum(pattern3)
			suma4 = np.sum(pattern4)
			suma5 = np.sum(pattern5)
			suma6 = np.sum(pattern6)

			new_data1[j].append(suma1)
			new_data2[j].append(suma2)
			new_data3[j].append(suma3)
			new_data4[j].append(suma4)
			new_data5[j].append(suma5)
			new_data6[j].append(suma6)


	matrix = np.zeros((len(data1), len(data1)))
	for i in range(len(new_data1)):
		for j in range(len(new_data1)):
			d1 = np.linalg.norm(np.array(new_data1[i]) - np.array(new_data1[j]))
			d2 = np.linalg.norm(np.array(new_data2[i]) - np.array(new_data2[j]))
			d3 = np.linalg.norm(np.array(new_data3[i]) - np.array(new_data3[j]))
			d4 = np.linalg.norm(np.array(new_data4[i]) - np.array(new_data4[j]))
			d5 = np.linalg.norm(np.array(new_data5[i]) - np.array(new_data5[j]))
			d6 = np.linalg.norm(np.array(new_data6[i]) - np.array(new_data6[j]))
			matrix[i,j] = d1+d2+d3+d4+d5+d6


	assignments = sch.fcluster(sch.linkage(matrix, method='ward'), 4, criterion='maxclust')
	assignments = np.array(assignments)-1 #Perque estiguin entre 0 i n_labels-1
	print(assignments)
	print(window_length, ':', check_assignments(assignments), ':', calinski_harabasz_score(matrix, assignments))
	








