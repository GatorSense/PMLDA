
clc;clear all;close all  %demo for pmlda

addpath support/lightspeed         % install lightspeed toolbox and add to your path
addpath support
%%% Generate Simulated Data

b=0.5;                        % Prior for the scaling factor s
N=[100,120,200]*2;            % generate three documents with 100, 120, 200 visual words,respectively
trueCluster{1}.mu=[8; 0];
trueCluster{1}.cov=eye(2);
trueCluster{2}.mu=[0; -8];
trueCluster{2}.cov=eye(2);
trueCluster{3}.mu=[-8; 0];
trueCluster{3}.cov=eye(2);
trueCluster{4}.mu=[0; 8];
trueCluster{4}.cov=eye(2);
topic=length(trueCluster);    % number of topics
alpha=ones(1,topic);          % Prior for topic proportion pi
for d=1:length(N)
    Data{d}=genSynData(N(d),alpha,b,trueCluster);
end


para.topic=topic;    % number of topics
para.alpha=alpha;    % Prior for topic proportion pi
para.exp=b;          % 1/para.exp = expectation of the scaling factor s
para.iter=500;       % number of iterations for Gibbs sampling
para.muStep=eye(2);  % variance of the proposal distribution for cluster means

%[samples,cluster]=pmlda_leastsquare(Data,para);
[samples,cluster]=pmlda_gibbs(Data,para);

% Visualize the result
figure();
subplot(2,3,1);title('doc 1: true Z & centers');hold on;
scatter(Data{1}.X(:,1),Data{1}.X(:,2),[],samples(1).zStar(:,1)); colormap(jet); 
plot(trueCluster{1}.mu(1),trueCluster{1}.mu(2),'+m','LineWidth',4,'MarkerSize',12);hold on
plot(trueCluster{2}.mu(1),trueCluster{2}.mu(2),'+m','LineWidth',4,'MarkerSize',12);hold on
plot(trueCluster{3}.mu(1),trueCluster{3}.mu(2),'+m','LineWidth',4,'MarkerSize',12);hold on
plot(trueCluster{4}.mu(1),trueCluster{4}.mu(2),'+m','LineWidth',4,'MarkerSize',12);hold on
subplot(2,3,2);title('doc 2: true Z & centers');hold on;
scatter(Data{2}.X(:,1),Data{2}.X(:,2),[],samples(2).zStar(:,1)); colormap(jet); 
subplot(2,3,3);title('doc 3: true Z & centers');hold on;
scatter(Data{3}.X(:,1),Data{3}.X(:,2),[],samples(3).zStar(:,1)); colormap(jet); 

subplot(2,3,4);title('doc 1: est Z & centers');hold on;
scatter(Data{1}.X(:,1),Data{1}.X(:,2),[],Data{1}.Z(:,1)); colormap(jet);
plot(cluster.mu{1}(1),cluster.mu{1}(2),'+k','LineWidth',4,'MarkerSize',12);hold on
plot(cluster.mu{2}(1),cluster.mu{2}(2),'+k','LineWidth',4,'MarkerSize',12);hold on
plot(cluster.mu{3}(1),cluster.mu{3}(2),'+k','LineWidth',4,'MarkerSize',12);hold on
plot(cluster.mu{4}(1),cluster.mu{4}(2),'+k','LineWidth',4,'MarkerSize',12);hold on
subplot(2,3,5);title('doc 2: est Z & centers');hold on;
scatter(Data{2}.X(:,1),Data{2}.X(:,2),[],Data{2}.Z(:,1)); colormap(jet); 
subplot(2,3,6);title('doc 3: est Z & centers');hold on;
scatter(Data{3}.X(:,1),Data{3}.X(:,2),[],Data{3}.Z(:,1)); colormap(jet); 



