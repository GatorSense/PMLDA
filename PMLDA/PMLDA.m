%% LDA-PM, by Chao Chen, 
% Affiliation: Tigersense Lab, University of Missouri-Columbia
% Copyrights 2016
% Topic proportion - pi ~ Dir(alpha)
% Scaling factor - s ~ EXP(lamda) = lamda*exp(-lamda*s)
% Membership (each point) - s ~ Dir(s*pi)
% Lightspeed toolbox should be installed.


function [cluster,samples] = PMLDA(HSID,parameters)
% Input: HSI data structure - HSID(n) represents the nth document
%                             each document should be d x N. d is the
%                             dimension and N is number of points in the
%                             document.
% Example: 
%         load 'HSID.mat'
%         [parameters] = PMLDAParameters();
%         [cluster,samples] = PMLDA(HSID);


addpath support
addpath support/lightspeed


%% Model Setup
if nargin < 2
[parameters] = PMLDAParameters();
end


%% load data. Only data points are necessary

D = length(HSID); % number of documents
for d = 1:D
    Data(d).X = HSID(d).Data';
end

dim = size(Data(1).X,2);

XX = [];% the whole dataset
for d = 1:D
    X = Data(d).X;% dataset
    N = size(X,1);% number of points
    
    % random initialization for membership
    z0 = rand(N,parameters.Ntopic);
    z0 = z0./repmat(sum(z0,2),1,parameters.Ntopic);% randomly initialize membership vector and normalize it sum-to-one
    samples(d).z = z0;
    
    
    % random initialization for topic proportion
    pi0 = rand(1,parameters.Ntopic);
    pi0 = pi0/sum(pi0);% randomly initialize topic proportion vector
    samples(d).pi = pi0;
    
    % random initialization for scaling factor
    s0 = rand;% radomly initialize scaling factor
    samples(d).s = s0;
    XX = [XX;X];% whole dataset
end

muX = mean(XX,1);
covX = cov(XX);

for q = 1:parameters.Ntopic
    cluster(q).mu = muX; % initialize cluster centers
    cluster(q).cov = 0.1*eye(dim,dim); % initialize cluster covariances
    newMu(q,:) = muX; % initialize cluster centers
end

Picandidate = dirichlet_sample(parameters.alpha,parameters.Iteration); % topic proportion candidates
Scandidate = exprnd(parameters.s,parameters.Iteration,1); % scaling factor candidates

num = parameters.Ntopic;
%% Main - PMLDA iterations
for iter = 2:parameters.Iteration
    iter
    
    allZ = [];
    for d = 1:D
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%Sample the membership value
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        s = samples(d).s;
        pi = samples(d).pi;
        X = Data(d).X;
        N = size(X,1);% number of points in each document
        for q = 1:parameters.Ntopic
            Nat{q,1} = inv(cluster(q).cov);
            Nat{q,2} = cluster(q).mu*Nat{q,1};
        end
        
        for n = 1:N
            newZn = dirichlet_sample(ones(1,parameters.Ntopic)*1,1);
            [newP] = compP(newZn,Nat,num,X,dim,n,s,pi);
            oldZn = samples(d).z(n,:);
            [oldP] = compP(oldZn,Nat,num,X,dim,n,s,pi);
            az = min(1,exp(newP - oldP));
            u = rand;
            if u < az
                samples(d).z(n,:) = newZn;
                ez = newP;
            else
                ez = oldP;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%Sample the topic proportion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        z = samples(d).z; %update z
        
        newPi = Picandidate(iter,:);
        newP = sum((parameters.alpha - 1).*log(newPi)) + N*gammaln(s) - N*sum(gammaln(s*newPi)) + sum(log(z)*(s*newPi' - 1));
        oldPi = samples(d).pi;
        oldP = sum((parameters.alpha - 1).*log(oldPi)) + N*gammaln(s) - N*sum(gammaln(s*oldPi)) + sum(log(z)*(s*oldPi' - 1));
        api = min(1,exp(newP - oldP));
        u = rand;
        if u < api
            samples(d).pi = newPi;
            ePi = newP;
        else
            ePi = oldP;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%Sample the scaling factor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pi = samples(d).pi;
        newS = Scandidate(iter);
        newP = -1/parameters.s*newS + N*gammaln(newS) - N*sum(gammaln(newS*pi)) + sum(log(z)*(newS*pi' - 1));
        oldS = samples(d).s;
        oldP = -1/parameters.s*oldS + N*gammaln(oldS) - N*sum(gammaln(oldS*pi)) + sum(log(z)*(oldS*pi' - 1));
        as = min(1,exp(newP - oldP));
        u = rand;
        if u < as
            samples(d).s = newS;
            es = newP;
        else
            es = oldP;
        end
        allZ = [allZ;z];
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Sample the clusters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for wave = 1:dim
        for q = 1:num
            Mu(q,:) = cluster(q).mu;
            Sig{q} = cluster(q).cov;
        end
        
        %%% update each cluster
        
        for q = 1:num
            newMu=cluster(q).mu;    
            newMu(wave)=mvnrnd(cluster(q).mu(wave),covX(wave,wave),1);
            
            oldtemp=[];
            newtemp=[];
            for r = 1:num
                oldtemp = horzcat(oldtemp,logmvnpdf(XX,Mu(r,:),Sig{r})');
                if r == q
                    newtemp = horzcat(newtemp,logmvnpdf(XX,newMu,Sig{r})');
                else
                    newtemp = horzcat(newtemp,logmvnpdf(XX,Mu(r,:),Sig{r})');
                end
            end
            
            newP = allZ.*newtemp;
            newP = sum(newP(:));
            oldP = allZ.*oldtemp;
            oldP = sum(oldP(:));
            
            aC = min(1,exp(newP - oldP));
            u = rand;
            if u < aC
                cluster(q).mu = newMu;
                eMu(q) = newP; 
                cStar.eng = newP;
            else
                eMu(q) = oldP;
            end
            
            Mu(q,:) = cluster(q).mu;
        end
                
    end
    
     %sample the covariance matrix    

     for r = 1:num
         newSig{r} = diag(ones(1,dim)*parameters.Sigcandidate(iter));
     end
     
     oldtemp = [];
     newtemp = [];
     for r = 1:num
         oldtemp = horzcat(oldtemp,logmvnpdf(XX,Mu(r,:),Sig{r})');
         newtemp = horzcat(newtemp,logmvnpdf(XX,Mu(r,:),newSig{r})');
     end
     newP = allZ.*newtemp;
     newP = sum(newP(:));
     oldP = allZ.*oldtemp;
     oldP = sum(oldP(:));
     
     aC = min(1,exp(newP-oldP));
     u = rand;
     if u < aC
         for r = 1:num
         cluster(r).cov = newSig{r};
         end
         cStar.eng = newP;
     end
     cStar.eng
     figure(10000)
for q = 1:num
    plot(1:dim,cluster(q).mu);title('estimated cluster centers');hold on;
end
pause(eps)
end

figure
for q = 1:num
    plot(1:dim,cluster(q).mu);title('estimated cluster centers');hold on;
end
hold off
end