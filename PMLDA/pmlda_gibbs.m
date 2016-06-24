function [samples,cluster]=pmlda_gibbs(Data,para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this function, both inference and cluster estimation are done using Gibbs sampler
% function "pmlda_leastsquare" updates the cluster using least square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input: Data - 1*D cells. Each cell represents a document.
%       Data{i}.X  - saves the visual words for document i, which is a
%                    N-by-p matrix, N is the number of visual word in document i, and p
%                    is the feature dimension.
%                    Note that users can add additional fields to Data{i}. This code only use
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
%
% Written by Chao Chen, University of Missouri-Columbia, 06/13/2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XX=[];
Length=[];
for i=1:length(Data)
    XX=[XX; Data{i}.X];
    Length=[Length size(Data{i}.X,1)];
end


dim=size(XX,2); % feature length
OFFSET=cumsum(Length)-Length; % saves the number of visual words in each document

%%%%%Parameter Settings%%
D=length(Data);
K=para.topic;
alpha=para.alpha;
b=para.exp;
Iteration=para.iter;
sampleStep=para.muStep;
%%%%%%Initialize the clusters using kmeans%%

%[c0,~]=vl_kmeans(XX',K);
[~, c0]=kmeans(XX,K);
c0=c0';
allMu=[];
allCov=[];
for k=1:K
    cluster.mu{k}=c0(:,k);
    cluster.cov{k}=eye(dim);
    cluster.nat1{k}=inv(cluster.cov{k});
    cluster.nat2{k}=cluster.nat1{k}*cluster.mu{k};
    allMu=[allMu; cluster.mu{k}(:)'];
    allCov=[allCov; cluster.cov{k}(:)'];
end

%%%%Initialize membership z, topic proportion pi, scaling factor s

for d=1:D
    N=Length(d);
    z0=rand(N,K);
    z0=z0./repmat(sum(z0,2),1,K);
    samples(d).z=z0;
    samples(d).zStar=z0;

    pi0=rand(1,K);
    pi0=pi0/sum(pi0);
    samples(d).pi=pi0;
    samples(d).piStar=pi0;

    s0=rand;
    samples(d).s=s0;
    samples(d).sStar=s0;

    energyZ=zeros(N,1);
    offset=OFFSET(d);
    for n=1:N
        Nat1=zeros(dim,dim);
        Nat2=zeros(dim,1);
        for k=1:K
           Nat1=Nat1+z0(n,k)*cluster.nat1{k};
           Nat2=Nat2+z0(n,k)*cluster.nat2{k};
        end
        covr=inv(Nat1);
        mu=covr*Nat2;
        eZ=logmvnpdf1(XX(offset+n,:),mu',covr);  
        eZ=sum(eZ(:));
        eZ=eZ+log(z0(n,:))*(s0*pi0'-1);
        energyZ(n,1)=eZ;
    end

    samples(d).zStarEng=energyZ;
    samples(d).sStarEng=-b*s0+N*gammaln(s0)-N*sum(gammaln(s0*pi0))+sum(log(z0)*(s0*pi0'-1));
    samples(d).piStarEng=sum((alpha-1).*log(pi0))+N*gammaln(s0)-N*sum(gammaln(s0*pi0))+sum(log(z0)*(s0*pi0'-1));

end



for iter=2:Iteration
    
    if(mod(iter,10)==0)
    display(['Iter: ' num2str(iter)] );
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%Sample the membership value
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    allZ=[]; % concatenation of memberships for all documents
    for d=1:D
        s=samples(d).s;
        pi=samples(d).pi(1,:);
        offset=OFFSET(d);
        N=Length(d);
        
      for n=1:N
          newZn=dirichlet_sample(ones(1,K),1);
          Nat1=zeros(dim,dim);
          Nat2=zeros(dim,1);
         for k=1:K
             Nat1=Nat1+newZn(k)*cluster.nat1{k};
             Nat2=Nat2+newZn(k)*cluster.nat2{k};
         end
         covr=inv(Nat1);
         mu=covr*Nat2;
         newP=logmvnpdf1(XX(offset+n,:),mu',covr);        
         newP=sum(newP(:));
         newP=newP+log(newZn)*(s*pi'-1);
       
         oldZn=samples(d).z(n,:);
         Nat1=zeros(dim,dim);
         Nat2=zeros(dim,1);
         for k=1:K
             Nat1=Nat1+oldZn(k)*cluster.nat1{k};
             Nat2=Nat2+oldZn(k)*cluster.nat2{k};
         end
         covr=inv(Nat1);
         mu=covr*Nat2;
         oldP=logmvnpdf1(XX(offset+n,:),mu',covr);        
         oldP=sum(oldP(:));
         oldP=oldP+log(oldZn)*(s*pi'-1);
       
         az=min(1,exp(newP-oldP));
         u=rand;
         if u<az
            samples(d).z(n,:)=newZn;
            ez=newP;
         else
            ez=oldP;
         end
      
         if ez>samples(d).zStarEng(n)  
            samples(d).zStar(n,:)=samples(d).z(n,:);
            samples(d).zStarEng(n)=ez;            
        end
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%Sample the topic proportion
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    z=samples(d).z; 
    newPi=dirichlet_sample(ones(1,K)*1,1);
    newP=sum((alpha-1).*log(newPi))+N*gammaln(s)-N*sum(gammaln(s*newPi))+sum(log(z)*(s*newPi'-1));
    oldPi=samples(d).pi(1,:);
    oldP=sum((alpha-1).*log(oldPi))+N*gammaln(s)-N*sum(gammaln(s*oldPi))+sum(log(z)*(s*oldPi'-1));
    api=min(1,exp(newP-oldP));
    u=rand;
    if u<api
       samples(d).pi=newPi;
       ePi=newP;
    else
       ePi=oldP;
    end
    if ePi>samples(d).piStarEng
       samples(d).piStar=samples(d).pi;
       samples(d).pStarEng=ePi;
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%Sample the scaling factor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pi=samples(d).pi;
    newS=exprnd(1/b,1);
    newP=-b*newS+N*gammaln(newS)-N*sum(gammaln(newS*pi))+sum(log(z)*(newS*pi'-1));
    oldS=samples(d).s;
    oldP=-b*oldS+N*gammaln(oldS)-N*sum(gammaln(oldS*pi))+sum(log(z)*(oldS*pi'-1));
    as=min(1,exp(newP-oldP));
    u=rand;
    if u<as
       samples(d).s=newS;
       es=newP;
    else
      es=oldP;
    end
  
    if es>samples(d).sStarEng
       samples(d).sStar=newS;
       samples(d).sStarEng=es;
    end
    allZ=[allZ;z];
 end
   

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% update the clusters
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for k=1:K
     newMu=mvnrnd(cluster.mu{k}',sampleStep,1);
     oldP=newProbability1(XX,allZ,allMu,allCov,K);
     temp=allMu(k,:);
     allMu(k,:)=newMu(:)';
     newP=newProbability1(XX,allZ,allMu,allCov,K);
     aC=min(1,exp(newP-oldP));
     u=rand;
     if u>aC
        allMu(k,:)=temp;
     end
 end
 for k=1:K
    cluster.mu{k}=allMu(k,:)';
    cluster.nat2{k}=cluster.nat1{k}*cluster.mu{k};   
 end
 
end





