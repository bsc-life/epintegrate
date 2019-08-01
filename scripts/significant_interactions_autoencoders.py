import pickle
import numpy as np
import random
from keras.layers import Input, Dense
from keras.models import Model
import pandas as pd
import os
import csv
import pickle

num=20
num_samples = 100
elim_number = 100
square_length = 50
elim_pos = []
importance_matrix = np.zeros((2490, 2490))

x = []

for exp in ['CLL_12_', 'CLL_110_', 'CLL_1525_', 'NBC_', 'MBC_', 'GCBC_']:
        print(exp)
        for i in range(1,num+1):
                print(i)
                path = 'hic_matrices/' + exp + str(i).zfill(2) + '/05_sub-matrices/'
                for file in os.listdir(path):
                        if file.endswith('.mat'):
                                fh = open(path + file)
                                next(fh)
                                next(fh)
                                interac_matrix = []
                                for line in fh:
                                        v = line.split()
                                        y = v[3:len(v)]
                                        y = map(float,y)
                                        interac_matrix.append(y)
                x.append(interac_matrix)


for i in range(num_samples):
	print(i)
	xx = list(x)
	xx = np.asarray(xx)

	for j in range(elim_number):
		row = random.randint(0, int(2490/square_length))
		col = random.randint(0, int(2490/square_length))
		elim_pos.append(tuple([row, col]))


	for j in range(len(xx)):
		for row, col in elim_pos:
			for k in range(square_length):
				for l in range(square_length):
					xx[j][row:row+k][row:row+l] = 0.0

	encoding_dim = 2490*2490 // 4

	for i in range(len(xx)):
		xx[i] = xx[i].reshape(2490*2490, 1)
	
	input_img = Input(shape=(2490*2490,))
	encoded = Dense(encoding_dim, activation='relu')(input_img)
	decoded = Dense(2490*2490, activation='sigmoid')(encoded)
	autoencoder = Model(input_img, decoded)
	#encoder = Model(input_img, encoded)

	#encoded_input = Input(shape=(encoding_dim,))
	#decoder_layer = autoencoder.layers[-1]
	#decoder = Model(encoded_input, decoder_layer(encoded_input))

	autoencoder.compile(optimizer='adadelta', loss='binary_crossentropy')

	history = autoencoder.fit(xx, xx, epochs=50, shuffle=True)

	loss = history.history['loss'][len(history.history['loss'])-1]

	for row, col in elim_pos:
		importance_matrix[row][col] += loss


np.savetxt('importance_matrix', importance_matrix)





	




