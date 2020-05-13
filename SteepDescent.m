%mnguye62 MATLAB R2019b
function [data_xyz, time, U, tol] = SteepestDescent(data_xyz,seq,hp_all,k,le,dt,m,temp,zeta)

time = zeros(1,1E6);
U = zeros(1,1E6);
count = 1;
F0= zeros(length(data_xyz),3);
[ U(1,2) , F1, ~] = nppm(data_xyz,seq,hp_all,k,le,temp,m,zeta,dt);
tol =zeros(1,length(1E6));


while sum(vecnorm((F1-F0),2),'all') > 1E-6
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
   [U(1,count), F1, ~] = nppm(data_xyz,seq,hp_all,k,le,temp,m,zeta,dt);
   tol(1,count) = sum(vecnorm((F1-F0),2),'all');

end

U = U(1,1:count);
time = time(1,1:count);

end