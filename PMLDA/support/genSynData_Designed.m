function data = genSynData_Designed(alpha,b,cluster,s)

% inputs:
% n - number of data points to be generated
% alpha - 1xk vector of dirichlet hypers
% b - scalar hyper for exponential distribution
% cluster - cluster parameters
% Chao 



k=length(alpha);
pi=dirichlet_sample(alpha); %topic proportion
%s=gamrnd(repmat(pi,n,1),1); %each data point has a different scale
m=1;
S=[];
Z=[];
XX=[];

if k==2
   D1=[ones(800,1)  zeros(800,1)];
   D2=[zeros(800,1) ones(800,1)];
   D3=rand([1600,1]);
   D3=[1-D3 D3];
   D=[D1;D2;D3];
end

if k==3
   D1=[ones(200,1)  zeros(200,1) zeros(200, 1)];
   D2=[zeros(200,1) ones(200,1) zeros(200, 1)];
   D3=[zeros(200,1) zeros(200,1) ones(200,1)];
   D12=rand([800,1]);
   D12=[1-D12 D12 zeros(800, 1)];
   D23=rand([800,1]);
   D23=[zeros(800, 1) 1-D23 D23];
   D13=rand([800,1]);
   D13=[1-D13  zeros(800, 1) D13];
   D=[D1;D2;D3;D12;D23;D13];    
end

%z=[0.2 0.8];
for ii=1:size(D,1)
    %z=dirichlet_sample(s*pi,1);
    z=D(ii,:);
%     [r c]=find(z==max(z));
%     z=0*z;
%     z(r,c)=1;
%     z
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
       S=[S ;s];
        
end

    data.X=XX;
    data.Z=Z;
    data.Scale=S;
    data.Pi=pi;
    data.cluster=cluster;
    