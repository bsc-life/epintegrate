import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
from sklearn.cluster import KMeans

df = pd.read_csv('more_data2.csv', sep=' ', header=None)
x = df.values

num_classes = 7
genome_length = 2490

def calinski_harabasz_score(X, labels):
	n_samples, _ = X.shape
	n_labels = len(set(labels))

	extra_disp, intra_disp = 0., 0.
	mean = np.mean(X, axis=0)
	
	for k in range(n_lables):
		cluster_k = X[labels == k]
		mean_k = np.mean(cluster_k, axis=0)
		extra_disp += len(cluster_k)*np.sum((mean_k - mean) ** 2)
		intra_disp += np.sum((cluster_k - mean_k) ** 2)

	return (1. if intra_disp == 0. else extra_disp * (n_samples - n_labels) / (intra * (n_labels - 1.)))

def get_importance_vector(i, j):
	xlist = list(x)
	data = xlist[20*i:20*(i+1)] + xlist[20*j:20*(j+1)]
	data = np.array(data)
	importance_vector = np.zeros(genome_length)
	occurrence_vector = np.zeros(genome_length)
	num_samples = 100000
	remaining_percentage = 0.01
	remaining_number = int(genome_length*remaining_percentage)

	for i in range(num_samples):
		if (i % 10000 == 0): print(i)
		remaining_positions = []
		xx = []

		for i in range(remaining_number):
			index = random.randint(0, len(data[0])-1)
			remaining_positions.append(index)

		for i in range(len(remaining_positions)):
			occurrence_vector[remaining_positions[i]] += 1

		for i in range(len(data)):
			v = [data[i][j] for j in remaining_positions]
			xx.append(v)

		xx = np.array(xx)
	
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

	return importance_vector	
	
x_mean = [np.mean(x[20*i:20*(i+1), axis=0) for i in range(num_classes)]
x_mean = np.asarray(x_mean)

def difference_to_others(vector2):
	aux = 0
	for i in [num for num in range(num_classes) if (num != 1 and num != 2)]:
		aux += np.linalg.norm(vector2 - x_mean[i])
	return aux

num = 20
order_matrix = np.zeros((num_classes, num_classes, num))

for k in range(num_classes):
	for l in range(num_classes):
		imp_vector = get_importance_vector(k, l)
		imp_regions = imp_vector.argsort()[-num:][::-1]
		used_regions = np.zeros(num)
		vector = xx[k]
		order = []
		for i in range(num):
			index = -1
			best = -np.inf
			for j in [d for d in range(num) if used_regions[d] == 0]:
				vector2 = vector
				vector2[imp_regions[j]] = xx[l][imp_regions[j]]
				aux = difference_to_others(vector2)
				if (aux > best):
					best = aux
					index = j
			used_regions[index] = 1
			vector[imp_regions[j]] = x_mean[l][imp_regions[j]]
			order.append(index)

		order_matrix[k][l] = order
	


























 
