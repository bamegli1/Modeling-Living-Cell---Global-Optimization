function[V]=graddescent(pos)

pmatrix = [1,3,9,12,14,16,17,18,21,22,23,24,25,26];
sigma = 1;
epsilon = 1;

[fX,fY,fZ,V]=energyFunction(sigma,epsilon,pos,pmatrix);
force = [fX.',fY.',fZ.'];
dt = 0.003;

while (true)
    
    postemp = pos + force*dt;
    [fX,fY,fZ,tempV]=energyFunction(sigma,epsilon,postemp,pmatrix);
    force = [fX.',fY.',fZ.'];
    
    pos = postemp;
    
    if abs(tempV-V) < 1e-6
        break
    end   
    V = tempV;

end

end