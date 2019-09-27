# Variants

### scriptFilterGWAS.R
Simple script to filter GWAS catalog variants. It searches for variants reported more than once (different studies, etc.) and only keeps the one with the lowest significance. The filtered GWAS catalog variants dataset is saved into a new file.

### scriptJoinVariants.R
Reads files containing variants associated to the disease from GWAS catalog, DisGeNET and ClinVar, and joins them in a unique file variants.data.csv. Mainly, GWAS catalog is considered; DisGeNET and ClinVar are used to extend the variants dataset with variants not present in the GWAS catalog.

# Loops

### findLoopsHomer.sh
Uses Homer Software to search significant chromatin contacts from a TADbit HiC-BAM file. It includes the step of converting the TADbit HiC-BAM file to Hi-C Summary format, for which having the BAMtoHiCsummary.py script in the same working directory is required. It takes the following arguments: input BAM file, output directory where to store the results, and Homer options (number of CPUs, resolution, window, genome assembly, BED file with “bad regions”).

### BAMtoHiCsummary.py
Takes a TADbit HiC-BAM file as input and converts it to the specialized Hi-C Summary Format for Homer (http://homer.ucsd.edu/homer/interactions/HiCtagDirectory.html). It takes the following arguments: input TADbit HiC-BAM file, output directory where to store the Hi-C Summary file and number of CPUs to use to read the BAM file.

### findLoopsHiCCUPS.sh
Uses HiCCUPS algorithm from Juicer to search significant chromatin contacts from a TADbit HiC-BAM file. It includes the step of converting the TADbit HiC-BAM file to 4DN DCIC format, for which having the BAMto4DNDCIC.py script in the same working directory is required. It takes the following arguments: input BAM file, header for 4DN DCIC format, output directory to store the results, genome (e.g. hg38) and HiCCUPS options (maximum size of the submatrix within the chromosome, chromosome(s) on which HiCCUPS will be run, resolution(s) for which HiCCUPS will be run, number of CPUs).

Header 4DN DCIC format example:  
\#\# pairs format v1.0  
\#columns: readID chr1 position1 chr2 position2 strand1 strand2  
\#chromsize: 1 248956422  
\#chromsize: 2 242193529  
\# ...  

### BAMto4DNDCIC.py
Takes a TADbit HiC-BAM file as input and converts it to the 4DN DCIC format for HiCCUPS (https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md). It takes the following arguments: input TADbit HiC-BAM file, input 4dn format header with pairs format, columns and chromsize, output directory where to store the 4DN file and number of CPUs to use to read the BAM file.

### Matrix_loops_comparison.ipynb
Compares the loops from two different BED files using a metric of distance and plots the comparison as a matrix. It also plots histograms of most similar loops (<3kb distance, <15kb distance).

### scriptJoinAllLoops.R
Simple script to read the loops found by Homer and HiCCUPS tools from their corresponding output directory and join them in a single file AllLoops.csv.

# Significant regions

### Check_length_ev.ipynb
Reads the eigenvectors obtained from TADbit and checks that the lengths correspond to those of the chromosomes from the Human Genome Assembly GRCh38 (https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38) in the specified resolution (in the default script 100kb).

### Join_genome_v2.ipynb
Reads the eigenvectors obtained from TADbit for all cell lines, all replicates and all chromosomes and writes them in a single file for the complete genome. Applies sigmoid transformation and ensures that all the eigenvectors for the complete genome have the same length by reaching the chromosome lengths in the Human Genome Assembly GRCh38. The output file format is a pandas dataframe with the following columns: cell line, replicate, chromosome, position and value.

### sigreg_pairwise_genome.py
Adaptation of sigreg_pairwise.py for the complete genome using output file from Join_genome_v2.ipynb.

### Join_genome_v3.ipynb
(Similar to Join_genome_v2.ipynb, changes the format of the output file).
Reads the eigenvectors obtained from TADbit for all cell lines, all replicates and all chromosomes and writes them in a single file for the complete genome. Applies sigmoid transformation and ensures that all the eigenvectors for the complete genome have the same length by reaching the chromosome lengths in the Human Genome Assembly GRCh38. The output file contains one line per cell line and replicate, containing the eigenvectors of all chromosomes, and being all of the same length.

### sigreg_pairwise_genome_original.py
Adaptation of sigreg_pairwise.py for the complete genome using output file from Join_genome_v3.ipynb.

# Integration

### scriptLoopsVariants.R
Reads the output files with significant loops  from Homer and HiCCUPS tools and finds which ones have variants falling in them. Variants file is a dataframe with at least chromosome ID and chromosomic position. In the default script, it is obtained from scriptJoinVariants.R. The intersecting loops with variants are saved in a table, LoopsVariants.csv.

### scriptPlotExtension.R
Plots different relationships between unique SNPs in loops and loops depending on KB extended in loops regions: 

* Number of loops vs. KB extended in loops regions
* Number of unique SNPs vs. KB extended in loops regions
* Number of unique SNPs from GWAS catalog database vs. KB extended in loops regions
* Number of unique SNPs from GWAS catalog databse with significant p-value ( <= 8*10**(-8)) vs. KB extended in loops regions

### scriptAllGWAS.R
Randomization using all available variants in the GWAS catalog, taking 244 (or number of variants associated to the disease) at random and counting how many loops have variants in their contact regions. 

### scriptComplementaryAnnotation.R
Reads table of intersecting loops with variants and annotates both the region of the loop with the variant and the region of the loop complementary to the variant. Annotates the regions without extending their ends (0 KB extended) and using the annotatr Bioconductor package. Writes a table with the loops regions where the variants fall and another table with the complementary regions. Also writes a table with the loops regions with the variants annotated and a table with the complementary regions annotated.

### scriptGlobalPlot_v2.R

### scriptCCscore.R

### scriptPlotProfiles.R
