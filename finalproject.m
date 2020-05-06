
function[]=finalproject()

data = importdata('onesequence_-22.79.dat',' ',1);
%P at 1,3,9,12,14,16,17,18,21,22,23,24,25,26
numid = [1;0;1;0;0;0;0;0;1;0;0;1;0;1;0;1;1;1;0;0;1;1;1;1;1;1;0];
pos = data.data;
Lambda=0.9999;
T = 100;
count = 0;

while T>0.001
    
    count = count + 1;
    [V,pos] = MCMC(3,T,pos);
    T = T*Lambda;
    
end
count
V
finalenergy = graddescent(pos)

%1 = P (polar), 0 = H (hydrophobic)
numid = [1;0;1;0;0;0;0;0;1;0;0;1;0;1;0;1;1;1;0;0;1;1;1;1;1;1;0]

writepdb(pos,numid,1,'27protein')
%step 1 
%file name string 


end

