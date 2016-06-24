Prerequisite:
Function newProbability1 and logmvpdf1 are implemented in C++, which require installation of the Eigen (a matrix operation) libary and GSL (GNU Scientific Library). These two libaries can be found in the following links:
http://eigen.tuxfamily.org/index.php?title=Main_Page
https://www.gnu.org/software/gsl/



The PM-LDA algorithm is realized in two ways.

1. [samples,cluster]=pmlda_gibbs(Data,para)

This function conducts both inference stage and cluster estimation stage using Gibbs sampler.

2. [samples,cluster]=pmlda_leastsquare(Data,para) 

This function conducts inference stage and cluster estimation stage using Gibbs sampler and least square method, respectively.

For the above two functions,

%Input: Data - 1*D cells. Each cell represents a document.
%       Data{i}.X  - saves the visual words for document i, which is a
%                    N-by-p matrix, N is the number of visual word in document i, and p
%                    is the feature dimension.
%                    Note that users can add additional fields to Data{i}. This function only uses
%                    Data{i}.X
%       para.topic - scalar, number of topics
%       para.alpha - 1-by-para.topic vector, parameter of Dirichlet distribution (prior of topic proportion)
%       para.exp   - scalar, parameter of exponential distribution (prior of scaling factor)
%       para.iter  - scalar, number of iterations for Gibbs sampling
%Output:
%      samples           - save the estimated memberships, topic proportion, scaling factor
%      samples(d).sStar  - estimated scaling factor for document d
%      samples(d).piStar - estimated topic proportion for document d
%      samples(d).zStar  - estimated membership for document d
%      cluster           - save the estimated clusters
%      cluster.mu{k}     - estimated mean for cluster k
%      cluster.cov{k}    - estimated covariance matrix for cluster k


Simulated data (for one document) can be generated using 

function data = genSynData(n,alpha,b,cluster)

% Inputs:
%        n       - number of data points (words, visual words) to be generated 
%        alpha   - 1xk vector of dirichlet hypers
%        b       - scalar hyper for exponential distribution
%        cluster - cluster parameters
%        cluster{k}.mu - mean of the kth topic
%        cluster{k}.cov - covariance matrix of the kth topic
% Outputs:
%        data    - a struct
%        data.X  - words, n-by-p matrix. p is the feature dimension
%        data.Z  - memberships, n-by-k matrix. k is the number of topics
%        data.Scale - scaling factor 
%        data.Pi - topic proportion
%        data.cluster - topics (clusters)



Please run demo.m to see how to run PM-LDA and visualize the estimated result. If you have any questions, please contact:

Chao Chen
Electrical and Computer Engineering
University of Missouri
Columbia, MO 65211
ccwwf@mail.missouri.edu

This code uses the Lightspeed Matlab toolbox (by Tom Minka), which can be downloaded from http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/
