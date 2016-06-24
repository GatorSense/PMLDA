# PMLDA
Partial Membership Latent Dirichlet Allocation

****************************************************************

NOTE: If PMLDA used in any publication or presentation, the following reference must be cited:  
C. Chen, A. Zare, and J. T. Cobb, "Partial Membership Latent Dirichlet Allocation for Image Segmentation," Under Review. 
available: http://arxiv.org/abs/1511.02821

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

Input: Data - 1*D cells. Each cell represents a document.  
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

 Chao Chen  
 Email Address:ccwwf@mail.missouri.edu  
 Alina Zare  
 Email Address:zarea@missouri.edu   
 University of Missouri, Department of Electrical and Computer Engineering  

****************************************************************

This product is Copyright (c) 2016 C. Chen, A. Zare  
 All rights reserved.  

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:  

   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.  
   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.  
   3. Neither the name of the University nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.  

 THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY OF MISSOURI AND
 CONTRIBUTORS ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,
 INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED.  IN NO EVENT SHALL THE UNIVERSITY OR CONTRIBUTORS
 BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES,
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE  


