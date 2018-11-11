# CSI-ANN

This package is developed with
OS: Linux 2.6.18-194.26.1.el5
Compiler: gcc version 4.1.2 20080704 (Red Hat 4.1.2-48).

Sample input files are provided in the Samples.zip file. 

The data is for human CD4+ T cell. References on the data source 
can be found at the end of this file. 

Extract sample files then follow the following steps for a try.

Step 1: Generate training and testing data
Compile: gcc -o s1 CSIANNPreprocess.c -lm
Run: ./s1

Input:
enhancers.txt  file containing information of known enhancers.
histones.txt   file containing histone modification data with a resolution of 200bp
genes.txt      gene annotation file

Output:
training.txt   file containing training data
testing.txt    file containing data for the prediction

Using the example files, it takes about 20 minutes to generate the training and testing 
data depending the speed of your computer. 


Step 2: Training the ANN model
Compile: make
Run: ./s2

Input: 
training.txt   training data generated in step 1    

Output:
Features.txt   FDA result 
partical_weights.txt    setting of the trained ANN model


Step 3: Do the prediction using the trained ANN model
Compile: gcc -o s3 CSIANNPredictionBatch.c -lm
Run: ./s3

Input:
testing.txt            testing data generated in step 1
Features.txt           FDA setting generated in step 2
partical_weights.txt   ANN setting generated in step 2

Output:
Prediction.txt         result of enhancer prediction


Please notice: All input files should be delimited by TAB.

Training enhancer markers (enhancers.txt) are derived from distal p300 binding sites
from the following reference:

Genome-wide mapping of HATs and HDACs reveals distinct functions in active and inactive genes.
Wang Z, Zang C, Cui K, Schones DE, Barski A, Peng W, Zhao K.
Cell. 2009 Sep 4;138(5):1019-31.

Histone modification data (.bed files) are obtained from the following reference:

Combinatorial patterns of histone acetylations and methylations in the human genome.
Wang Z, Zang C, Rosenfeld JA, Schones DE, Barski A, Cuddapah S, Cui K, Roh TY, Peng W, Zhang MQ, Zhao K.
Nat Genet. 2008 Jul;40(7):897-903.

Please cite the following reference for the CSI-ANN algorithm:

Discover regulatory DNA elements using chromatin signatures and artificial neural network.
Firpi HA, Ucar D, Tan K.
Bioinformatics. 2010 Jul 1;26(13):1579-86.