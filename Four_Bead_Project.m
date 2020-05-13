%% 
data = importdata('onesequence_-22.79.dat',' ',1);
xyz = data.data;
% Sequence 0 = nonpolar, 1 = polar, 2 = - charge, 3 = + charge
% arbitrarily Selected which polar residues became charged or not
seq = [ 1   0   1   0   0   0   0   0   2   0   0   3   0   2   0   1   3   1   0   0   2   1   3   1   3   2   0 ];

%Charge values for polarizability of a neutral carbon reduced units, dipole moment of water, plus charge, minus
%charge
nppm_all = [0.3135, 0.227153, -1, 1, 1, 1, 2/3];
%%
[coords, time, U, tol] = SteepDescent(xyz,seq,nppm_all,20,1,0.004,1,2,0.05);
figure(1)
plot(time,U)
figure(2)
plot3(coords(:,1),coords(:,2),coords(:,3),'-bo','LineWidth',5)
hold on 
%plot3(xyz(:,1),xyz(:,2),xyz(:,3),'-ro')
hold off
%%
[coords, time, ke, U] = nppmLD(xyz,seq,nppm_all,20,1,0.003,10000,1,2,0.05);
[coords, time, ke, U] = nppmLD(coords,seq,nppm_all,20,1,0.003,10000,1,2,0.05);
figure(1)
plot(time,ke,time,U,time,U+ke)
legend('KE','U','Total Internal Energy')
xlabel('Time (s)','FontSize',20)
ylabel('Energy (Epsilon)','FontSize',20)
title('Langevin Dynamics of Minimalist HP Protein Model','FontSize',20)
ave_ke = mean(ke(1,length(ke)-1000:length(ke)))
ave_pe = mean(U(1,length(U)-1000:length(U)))
figure(2)
plot3(coords(:,1),coords(:,2),coords(:,3),'-bo')
hold on 
plot3(xyz(:,1),xyz(:,2),xyz(:,3),'-ro')
%plot3(coords1(:,1),coords1(:,2),coords1(:,3),'-go')
hold off
%%
data = importdata('onesequence_-22.79.dat',' ',1);
xyz = data.data;
seq = cell2mat(data.textdata(2:length(data.textdata),:));
seq = seq(:,1);
hp_all = [1,2/3,2/3];
%% 4.a.ii
data = importdata('onesequence_-22.79.dat',' ',1);
xyz = data.data;
seq = cell2mat(data.textdata(2:length(data.textdata),:));
hp_all = [1,2/3,2/3];
[coords, time, ke, U] = hpLD(xyz,seq,hp_all,1,20,1,0.003,10000,1,2,0.05);
[coords1, time, ke, U] = hpLD(coords,seq,hp_all,1,20,1,0.003,10000,1,2,0.05);
figure(20)
plot(time,ke,time,U,time,U+ke)
legend('KE','U','Total Internal Energy')
xlabel('Time (s)','FontSize',20)
ylabel('Energy (Epsilon)','FontSize',20)
title('Langevin Dynamics of Minimalist HP Protein Model','FontSize',20)
%%
figure(10)
plot3(xyz(:,1),xyz(:,2),xyz(:,3),'-ro')