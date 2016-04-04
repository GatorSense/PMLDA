%
function [parameters] = PMLDAParameters()

parameters.Ntopic = 2;% number of topics, i.e. number of clusters
parameters.alpha=10*ones(1,parameters.Ntopic);% alpha represents the parameter vector for dirichlet hypers.
parameters.s=2; % scaling factor
parameters.Iteration=10; % number of iterations
parameters.Sigcandidate = 0.4*ones(parameters.Iteration,1); % variance candidates. Note: the default value is dataset dependent instead of 0.4.