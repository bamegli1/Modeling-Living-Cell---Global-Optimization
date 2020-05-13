
function[]=finalproject()

data = importdata('onesequence_-22.79.dat',' ',1);
%P at 1,3,9,12,14,16,17,18,21,22,23,24,25,26
numid = [1;0;1;0;0;0;0;0;1;0;0;1;0;1;0;1;1;1;0;0;1;1;1;1;1;1;0];
pos = data.data;
Lambda=0.9999;
T = 100;
count = 0;
num_part = 27;

nppm_all = [0.3135, 0.227153, -1, 1, 1, 1, 2/3];
k1=20;
le=1;
dt=0.003;
m=1;
temp=2;
zeta=0.05;
seq = [ 1   0   1   0   0   0   0   0   2   0   0   3   0   2   0   1   3   1   0   0   2   1   3   1   3   2   0 ];

total_pos = zeros(500000, num_part, 3);
while T>0.001
    
    count = count + 1;
    [V,pos] = MCMC(3,T,pos,seq,nppm_all,k1,le,temp,m,zeta,dt);
    T = T*Lambda;
    if (mod(count,100) == 0)
        total_pos(count/100, :, :) = pos;
    end
end
count
V

[finalenergy, grad_total_pos, grad_count] = graddescent(pos,seq,nppm_all,k1,le,temp,m,zeta,dt);

finalenergy

%1 = P (polar), 0 = H (hydrophobic)
numid = [1;0;1;0;0;0;0;0;1;0;0;1;0;1;0;1;1;1;0;0;1;1;1;1;1;1;0];

numFrames = floor(count/100); %sped up 100x
total_pos = total_pos(1:numFrames, :, :);
total_pos = permute(total_pos, [2, 1, 3]); %reshape in correct order
total_pos = reshape(total_pos, num_part*numFrames, 3);


grad_total_pos = permute(grad_total_pos, [2, 1, 3]); %reshape in correct order
grad_total_pos = reshape(grad_total_pos, num_part*grad_count, 3);

gen_xyzfile(total_pos, 27, 'filename', 'mcmc_total_pos');
gen_xyzfile(grad_total_pos, 27, 'filename', 'grad_total_pos');
writepdb(pos,numid,1,'27protein');

%step 1 
%file name string 


end

