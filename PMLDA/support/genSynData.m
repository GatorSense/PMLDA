function data = genSynData(n,alpha,b,cluster,m)

% inputs:
% n - number of data points to be generated
% alpha - 1xk vector of dirichlet hypers
% b - scalar hyper for exponential distribution
% cluster - cluster parameters
% m - fuzzifier
% Chao 

k=length(alpha);
pi=dirichlet_sample(alpha); %topic proportion
s=exprnd(1/b,1); %scaling factor
Z=[];
XX=[];

for ii=1:n
    z=dirichlet_sample(s*pi);
    if(k==3)
       natpara1=z(1)^m*inv(cluster.cov1)+z(2)^m*inv(cluster.cov2)+z(3)^m*inv(cluster.cov3);
       natpara2=z(1)^m*inv(cluster.cov1)*cluster.mu1+z(2)^m*inv(cluster.cov2)*cluster.mu2+z(3)^m*inv(cluster.cov3)*cluster.mu3;
    else
       natpara1=z(1)^m*inv(cluster.cov1)+z(2)^m*inv(cluster.cov2);
       natpara2=z(1)^m*inv(cluster.cov1)*cluster.mu1+z(2)^m*inv(cluster.cov2)*cluster.mu2;
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
    