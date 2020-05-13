%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mnguye62 MATLAB R2019b
%Inputs:
%data_xyz: (:,3) matrix of residue coordinates
%seq: bead identities  0 = nonpolar, 1 = polar, 2 = - charge, 3 = + charge
%nppm_all: harge values for polarizability of a neutral carbon reduced units, dipole moment of water, plus charge, minus
%charge, sigma, episilon for hh, episilon for hp
%k1: force constant of bead to bead bonds
%le: equilibrium length of the bonds
%temp: temperature of the NVT simulation
%zeta:  constant for BBK algorithmn
%dt: time step
%steps: number of steps for LD
%m: mass of atom
%Outputs:
%time: vector of all time steps
%U: Vector of PE at each time step
%coords: Final Coordinates of the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coords,  time, ke, U] = nppmLD(data_xyz,seq,nppm_all,k,le,dt,steps,m,temp,zeta)
x = data_xyz(:,1);
y = data_xyz(:,2);
z = data_xyz(:,3);

% initialize velocities
v = sqrt(temp/m)*randn(size(data_xyz));
meanv = mean(v);
v = v - meanv;

t = 0;
time = zeros(1,steps);
ke = zeros(1,steps);
U = zeros(1,steps);
count = 0;

while count < steps
   count = count + 1;
   time(count) = t;
   t = t + dt; 
   ke(1,count) = (1/2)*m*sum(v.^2,'all');
   
   if count == 1
   [U(1,count), F1, Fr] = nppm([x,y,z],seq,nppm_all,k,le,temp,m,zeta,dt);
   F1 = F1 + Fr - zeta * m * v;
   end
   
   %Update the positions
   for i = 1:length(data_xyz)
        x(i,1) = x(i,1) + dt* v(i,1)+ ((dt^2)/(2*m))*(F1(i,1));
        y(i,1) = y(i,1) + dt* v(i,2)+ ((dt^2)/(2*m))*(F1(i,2));
        z(i,1) = z(i,1) + dt* v(i,3)+ ((dt^2)/(2*m))*(F1(i,3));
   end
   
   %Force Calculations update
   [U(1,count), F2, Fr2] = nppm([x,y,z],seq,nppm_all,k,le,temp,m,zeta,dt);

   %Updated Velocities
   for i = 1:length(data_xyz)
       v(i,1) = (1/(1+(dt*zeta)/2))*(v(i,1)+(dt/(2*m))*F1(i,1)+(dt/(2*m))*(F2(i,1)+Fr2(i,1)));
       v(i,2) = (1/(1+(dt*zeta)/2))*(v(i,2)+(dt/(2*m))*F1(i,2)+(dt/(2*m))*(F2(i,2)+Fr2(i,2)));
       v(i,3) = (1/(1+(dt*zeta)/2))*(v(i,3)+(dt/(2*m))*F1(i,3)+(dt/(2*m))*(F2(i,3)+Fr2(i,3)));
   end
   
   %New force = Old force
   F1 = F2 + Fr2 - zeta * m * v;
   
end

coords = [x,y,z];