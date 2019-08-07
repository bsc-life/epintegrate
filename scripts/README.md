## Content of the scripts

#### find_optimal_resol_histones.py

The data used in this script can be found at */gpfs/projects/bsc08/shared_projects/epideco/data/BINARIZE_files*.  
Histone Marks data from Naive B Cell, Plasma Cell, Germinal Center B-Cell and Class Switched Memory B Cell have been used.
Two or three samples per cell type are available.
This script finds the optimal resolution to be used so that the clusters are well formed.
This is done by checking the accuracy of the clusters formed and by optimizing a cluster quality metric: Calinski-Harabasz Index

#### significant_regions_histones.py
By taking the optimal resolution found in *find_optimal_resol_histones.py*, this script finds an importance vector associated to the clustering of all the data.  
First we take the original data from */gpfs/projects/bsc08/shared_projects/epideco/data/BINARIZE_files* 
and turn it to the optimal resolution.  
Then we successively subsample the data randomly and check the quality of the clustering.  
We eventually get a genome-wide vector in which each component represents the "importance" of this region in clustering the data.

#### sigreg_pairwise.py
This script uses the first eigenvector of the hic-matrix of chromatin interactions for different cell types.
We used CLL_1525, CLL_110, CLL_12, MBC, NBC, GCBC and PBC data.  
A similar algorithm to the one used for histone marks is implemented to get an importance vector.
Note that the importance vectors might have different lengths because different resolutions might have been used.

#### significant_regions_autoencoders.py
Here we use the entire matrix of chromatin interactions.
We randomly add noise (or anulate) regions of the matrix and use the distorted matrix as input of an autoencoder,
that tries to reconstruct the original matrix.  
At the end of training, we use the loss function the autoencoder has optimized to check how well the reconstruction went.
Eventually we get a matrix of reconstruction performance for each interaction.

#### esotheric.py
Script that implements a probably esotheric idea to find out the order of the changes in epigenetic marks during a transition between two cell types. It is implemented using the Hi-C data we previously described.  
The idea is to choose the sequence that, at any timestep, the series of changes done up to that point don't result in the epigenetic profile being similar to a cell type different to the two involved in the transition.  
The one that maximizes this difference to the other cell types is the one we choose.  
We do that for every possible transition between two cell types.
