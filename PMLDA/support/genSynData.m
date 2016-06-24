function data = genSynData(n,alpha,b,cluster)

% Inputs:
%        n       - number of data points to be generated
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
% Chao 

K=length(alpha);
pi=dirichlet_sample(alpha); %draw topic proportion
s=exprnd(1/b,1);            %draw scaling factor
dim=size(cluster{1}.mu,1);
Z=[];
XX=[];

for ii=1:n
    z=dirichlet_sample(s*pi); %draw membership
    natpara1=zeros(dim,dim);
    natpara2=zeros(dim,1);
    for k=1:K
       natpara1=natpara1+z(k)*inv(cluster{k}.cov);
       natpara2=natpara2+z(k)*inv(cluster{k}.cov)*cluster{k}.mu;
    end
    cov=inv(natpara1);
    mu=inv(natpara1)*natpara2; 
    xx=mvnrnd(mu',cov,1); 
    XX=[XX ;xx];
    Z=[Z ;z];       
end

data.X=XX;
data.Z=Z;
data.Scale=s;
data.Pi=pi;
data.cluster=cluster;
    