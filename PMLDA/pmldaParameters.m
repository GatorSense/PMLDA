function [para] = pmldaParameters()


para.topic=3;    % number of topics
para.alpha=ones(1,para.topic);    % Prior for topic proportion pi
para.exp=0.5;          % 1/para.exp = expectation of the scaling factor s
para.iter=500;       % number of iterations for Gibbs sampling
para.muStep=eye(3);  % variance of the proposal distribution for cluster means