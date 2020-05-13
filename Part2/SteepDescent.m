%mnguye62 MATLAB R2019b
function [data_xyz, time, U, tol, count,Xtrajectory,Ytrajectory,Ztrajectory] = SteepestDescent(data_xyz,seq,hp_all,sigma,k,le,dt,m,temp,zeta)

time = zeros(1,1E6);
U = zeros(1,1E6);
count = 1;
F0= zeros(length(data_xyz),3);
[ U(1,2) , F1, ~] = hpLJ(data_xyz,seq,hp_all,sigma,k,le,temp,m,zeta,dt);
tol =zeros(1,length(1E6));

Xtrajectory = zeros(27,2707);
Ytrajectory = zeros(27,2707);
Ztrajectory = zeros(27,2707);

b = vecnorm((F1-F0),2);

while sum(b(:)) > 1E-6
   Xtrajectory(:,count) = data_xyz(:,1);
   Ytrajectory(:,count) = data_xyz(:,2);
   Ztrajectory(:,count) = data_xyz(:,3);
  
   count = count + 1;
   time(count) = count; 
   
   %Comparing Potential
   if U(count) > U(count-1)
       a = 0.9; 
   else
       a = 1.1;
   end
  
   lam = dt * a;
  
   %Update the positions
   data_xyz = data_xyz +F1 * lam;
   F0=F1;
   
   %Gradient calculation
   [U(1,count), F1, ~] = hpLJ(data_xyz,seq,hp_all,sigma,k,le,temp,m,zeta,dt);
   b = vecnorm((F1-F0),2); tol(1,count) = sum(b(:));

   
end

Xtrajectory(:,count) = data_xyz(:,1);
Ytrajectory(:,count) = data_xyz(:,2);
Ztrajectory(:,count) = data_xyz(:,3);
   
U = U(1,1:count);
time = time(1,1:count);

end

