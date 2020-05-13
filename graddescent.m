function[V, positions, i]=graddescent(pos,seq,nppm_all,k1,le,temp,m,zeta,dt)

pmatrix = [1,3,9,12,14,16,17,18,21,22,23,24,25,26];
sigma = 1;
epsilon = 1;



[V, force, ~]=nppm(pos,seq,nppm_all,k1,le,temp,m,zeta,dt);
%dt = 0.003;
positions = zeros(100000, 27, 3);
i = 0;
while (true)
    i = i + 1;
    postemp = pos + force*dt;
    [tempV, force, ~]=nppm(postemp,seq,nppm_all,k1,le,temp,m,zeta,dt);
    %[fX,fY,fZ,tempV]=energyFunction(sigma,epsilon,postemp,pmatrix);
    %force = [fX.',fY.',fZ.'];
    
    pos = postemp;
    positions(i, :, :) = pos;
    
    if abs(tempV-V) < 1e-6
        break
    end   
    V = tempV;

end

positions = positions(1:i, :, :);
end