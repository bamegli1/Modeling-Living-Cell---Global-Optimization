% andrew margolis, 5/12/20
clc;clear;close all;

%% Initialize

data = importdata('onesequence_-22.79.dat',' ',1);
xyz = data.data;
seq = cell2mat(data.textdata(2:length(data.textdata),:));
seq = seq(:,1);
hp_all = [1,2/3,2/3];

%% Markov Chain with Steepest Descent

[coords, time, U, tol, count,Xtrajectory,Ytrajectory,Ztrajectory] = SteepDescent(xyz,seq,hp_all,1,20,1,0.004,1,2,0.05);
figure(3)
U(1,length(U))
time(1,length(time))
plot(time,U)
hold on
% yline(-66.6735,'r');
%-66.6735 was used because it was seen as the lowest energy value
hold off
legend('U')
xlabel('Number of Steps','FontSize',10)
ylabel('Energy (Epsilon)','FontSize',10)
title(' Minimalist HP Protein Model','FontSize',10)

%% Plotting 
figure(4)
plot3(coords(:,1),coords(:,2),coords(:,3),'-or','MarkerSize',7.5,'LineWidth',2.5)
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Most Stable Configuration of Bead Protein Structure')
a = seq == 'P';
% writepdb(coords,a,1,'4d')

%% Cluster
Kclusters = 25
[features,states,P,T50,means,feat_dist] = Clustering(Xtrajectory,Ytrajectory,Ztrajectory,Kclusters)
% spectral cluster
Mc = 5
[Ptransnew, peqnew, clust_reassignments, eigvalue, evnew, taunew, K, time_states ]=spectral_cluster(T50,P,10,Mc);

%% Visualization
stateNames = ['State 1','State 2','State 3','State 4'];
reduced_reassignments = zeros(Mc,2);reduced_reassignments(1,:) = [1,5];
counter = 2;
out = sortrows(clust_reassignments,2);
steps = length(Xtrajectory)
S = cell(Mc,1);

for i = 1:Mc
    OGstate = out(i,1);
    structure = means(:,OGstate);
    for j = 1:steps
        check = feat_dist(:,j);
        if structure == check
            Coordinates = [Xtrajectory(:,j),Ytrajectory(:,j),Ztrajectory(:,j)];
            S{i,1} = Coordinates;
        end
    end
end
