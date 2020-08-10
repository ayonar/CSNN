# Compressed Sensing for Network Networks

> Code for "A Compressed Sensing Framework for Efficient Dissection of Neural Circuits." Yonar and Lee et. al. (2019) [*Nature Methods* 16, pages126–133](https://www.nature.com/articles/s41592-018-0233-6)

Compressed Sensing based method for identifying key neurons/units in biological and artificial neural networks


## Identify Key Neurons
> This part includes the code for the inference of key neurons from behavioral data obtained from optogenetic screen as well as the code for inferring the key neurons in artificial neural network.

#### Sparse Solution

This script infers key neurons in controlling speed of locomotion in *C. elegans* neural network using behavioral data and the measurement matrix.
> SparseSolution.ipynb


#### Compressed Sensing on nonlinear Artificial Neural Network

Application of compressed sensing framework on a nonlinear artificial neural network trained to recognize handwritten digits to assess feasibility of the method.
> CompressedSensing_ANN.ipynb

*This script loads a trained network used in the paper. You can re-train a new model if you wish by uncommenting the corresponding cells.*


## Robustness Analysis
> This part includes the false negative/positive rate simulations and Arch efficiency simulations using corrupted measurement matrices.

#### Simulations

 - "FalseNegPosRate_Simulatuions.m" calculates the false negative false positive rates 
in identifying the essential interneurons. 
 - "Arch_Efficiency_Simulations.m"  obtains 1000 sparse solutions using 1000 corrupt measurement matrices to determine robustness
of the solution to variations in the levels of archeorhodopsin expression.

#### False Neg/Pos Rate Simulations

> This code will analyze three measurement matrix sizes by default. To choose just one, change line 7 of `run.sh.`

#### Arch Efficiency Simulations

> 1000 corrupt matrices were used to check robustness of the solution to errors in measurement matrix.


```

```

#### Dependencies
The code is tested in Python 2.7 and 3.6 and Matlab2017b, and on MacOs 10.13 and Ubuntu 16.04

miniconda 4.3 + the libraries below should be enough to run this part.
dill==0.2.6 \
jupyter==1.0.0 \
Keras==2.1.3 \
h5py==2.8.0 \
matplotlib==2.2.2 \
scikit-learn==0.19.2 \
seaborn==0.8 \
tensorflow==1.10.0