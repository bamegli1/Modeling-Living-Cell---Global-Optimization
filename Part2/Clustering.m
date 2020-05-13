function[features,states,P,T50,means,feat_dist] = Clustering(Xtrajectory,Ytrajectory,Ztrajectory,Kclusters)

% clc;clear;close all
%  
% load Zvstime.dat; 
% Z = Zvstime;
% load Yvstime.dat; 
% Y = Yvstime;
% load Xvstime.dat; 
% X = Xvstime;

%% Feature array
n = size(Xtrajectory,1);
structdists = zeros(nchoosek(n,2),2);

counter = 1 
for i = 1:(n-1)
    for j = (i+1):n
        structdists(counter,:) = [i,j];
        counter = counter + 1;
    end
end
        
Ndists = length(structdists);
Nstructures = size(Xtrajectory,2);

test_dist = zeros(Ndists,Nstructures);

for i = 1:Nstructures
    for j = 1:Ndists
        b1 = structdists(j,1); b2 = structdists(j,2);
        dX = Xtrajectory(b1,i)-Xtrajectory(b2,i);
        dY = Ytrajectory(b1,i)-Ytrajectory(b2,i);
        dZ = Ztrajectory(b1,i)-Ytrajectory(b2,i);
        test_dist(j,i) = sqrt(dX^2+dY^2+dZ^2);
    end
end

%% Determine which features to use
feat_var = zeros(Ndists,1);

for i = 1:Ndists
    variance = variancefun(test_dist(i,:));
    feat_var(i) = variance;
end
[out,idx] = sort(feat_var);

% choosing sorted feature variances #249 to #351
features = zeros((351-249),2);
counter = 1;
for i = 249:351
    features(counter,:) = structdists(idx(i),:);
    counter = counter + 1;
end

% structures =

Nfeatures = length(features);
Nstructures = size(Xtrajectory,2);

feat_dist = zeros(Nfeatures,Nstructures);

for i = 1:Nstructures
    for j = 1:Nfeatures
        b1 = features(j,1); b2 = features(j,2);
        dX = Xtrajectory(b1,i)-Xtrajectory(b2,i);
        dY = Ytrajectory(b1,i)-Ytrajectory(b2,i);
        dZ = Ztrajectory(b1,i)-Ytrajectory(b2,i);
        feat_dist(j,i) = sqrt(dX^2+dY^2+dZ^2);
    end
end

[states,means] = kmeans(feat_dist.',Kclusters);
means = means.'

% 6 x 6 transition matrix

stepsize = 10;
T50 = trans_mat(states,stepsize,Nstructures,Kclusters);

% 3c 
% - transitions between the same state occur most frequently
% - states 1,2,and 6 have the highest probability of transitioning
% - states 3,4 have the lowest probability of transitioning


%% Probability of being in each state

P = zeros(Kclusters,1);

for i = 1:Nstructures
    P(states(i)) = P(states(i)) + 1;
end

P = P./sum(P);

%% Transition Probabilities at 100 steps apart

stepsize = 100
T100 = trans_mat(states,stepsize,Nstructures,Kclusters);
