# PMLDA
Partial Membership Latent Dirichlet Allocation

****************************************************************

NOTE: If this code is used, cite it: Chao Chen, Alina Zare, & Timotius Andrean Patrick Lagaunne. (2019, April 12). GatorSense/PMLDA v1.0 (Version v1.0). Zenodo. http://doi.org/10.5281/zenodo.2638296 [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2638296.svg)](https://doi.org/10.5281/zenodo.2638296)

NOTE: If PMLDA used in any publication or presentation, the following reference must be cited:  
C. Chen, A. Zare, H. Trinh, G. Omotara, J. T. Cobb, and P. Lagaunne, “Partial Membership Latent Dirichlet Allocation,” IEEE Trans. Image Process., vol. 26, pp. 5590-5602, Dec. 2017. 
available: https://arxiv.org/abs/1612.08936

****************************************************************
PRIOR TO RUNNING PMLDA:  
+ PMLDA code uses the following toolboxes and libraries: Lightspeed Matlab toolbox, Eigen Library and GSL (GNU Scientific Library)  
    * Lightspeed Matlab toolbox (by Tom Minka), which can be downloaded from http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/
    * Eigen Library: http://eigen.tuxfamily.org/index.php?title=Main_Page
    * GSL Library: https://www.gnu.org/software/gsl/
+ The Lightspeed toolbox and the Eigen library must both be downloaded and placed in the support/ folder
+ GSL Library must be downloaded and compiled/installed 
+ After installing the libraries and toolboxes, the functions support/newProbability1.cpp and logmvnpdf1.cpp must be compiled with the following commands:
    * mex logmvnpdf1.cpp
    * mex newProbability1.cpp
+ Add paths to lightspeed, Eigen, and support   

****************************************************************

The command to run the PMLDA algorithm:   

[samples,cluster]=pmlda_gibbs(Data,para)

or 

[samples,cluster]=pmlda_leastsquare(Data,para) 

pmlda_gibbs conducts both inference stage and cluster estimation stage using Gibbs sampler.  

pmlda_leastsquare conducts the inference stage and the cluster estimation stage using a Gibbs sampler and a least squares approach, respectively.


For the above two functions,

Input:  
+ Data - 1*D cells. Each cell represents a document.  
      + Data{i}.X  - saves the visual words for document i, which is a
                    N-by-p matrix, N is the number of visual word in document i, and p
                    is the feature dimension.
                    Note that users can add additional fields to Data{i}. This function only uses
                    Data{i}.X  
+ para.topic - scalar, number of topics  
+ para.alpha - 1-by-para.topic vector, parameter of Dirichlet distribution (prior of topic proportion)  
+  para.exp   - scalar, parameter of exponential distribution (prior of scaling factor)  
+ para.iter  - scalar, number of iterations for Gibbs sampling  
       
Output:  
      + samples           - save the estimated memberships, topic proportion, scaling factor  
     +  samples(d).sStar  - estimated scaling factor for document d  
      + samples(d).piStar - estimated topic proportion for document d  
     +  samples(d).zStar  - estimated membership for document d  
      + cluster           - save the estimated clusters  
      + cluster.mu{k}     - estimated mean for cluster k  
      + cluster.cov{k}    - estimated covariance matrix for cluster k  

Simulated data (for one document) can be generated using   

function data = genSynData(n,alpha,b,cluster)

Inputs:  
      +  n       - number of data points (words, visual words) to be generated   
      +  alpha   - 1xk vector of dirichlet hypers  
      +  b       - scalar hyper for exponential distribution  
      +  cluster - cluster parameters  
      +  cluster{k}.mu - mean of the kth topic  
      +  cluster{k}.cov - covariance matrix of the kth topic    


Outputs:  
      +  data    - a struct  
      +  data.X  - words, n-by-p matrix. p is the feature dimension  
      +  data.Z  - memberships, n-by-k matrix. k is the number of topics  
      +  data.Scale - scaling factor   
      +  data.Pi - topic proportion  
      +  data.cluster - topics (clusters)  


Please run demo.m to see how to run PM-LDA and visualize the estimated result.  

****************************************************************  

Files explanation:  
Latest Revision: June 2016

+ README.md -  this file  
+ demo.m - demo script to run pmlda on synthetic data
+ LICENSE - license for the code
+ pmlda_gibbs.m - gibbs sampler implementation of pmlda
+ pmlda_leastsquare.m - gibbs/least_squares implementation of pmlda
+ /support - folder containing helper functions for pmlda implementation
+ /support/dirichlet_sample.m - function to generate a sample from dirichlet distribution
+ /support/logmvnpdf1.cpp 
+ /support/newProbability1.cpp 
+ /support/genSynData.m  


****************************************************************  

For any questions, please contact: 
 
 Alina Zare  
 Email Address:azare@ece.ufl.edu   
 University of Florida, Department of Electrical and Computer Engineering  

****************************************************************
