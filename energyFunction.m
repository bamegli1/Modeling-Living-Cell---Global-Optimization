%Brandon Ameglio (bamegli1)
%Inputs = sigma, epsilon, box lengh, and a data set.
%Outputs = potential energy of data. 
%Energy funciton for finding LJP of a dataset based on distances. 
%LD energy function of a particle

function[fX,fY,fZ,tempP]=energyFunction(sigma,epsilon,data,pmatrix)

k = 20;
le = 1;
epsilonnb = 2/3;

fX = zeros(1,length(data));
fY = zeros(1,length(data));
fZ = zeros(1,length(data));
tempP = 0;
for i=1:length(data)
    for j=i+1:length(data)
        point = [data(i,1)-data(j,1),data(i,2)-data(j,2),data(i,3)-data(j,3)];
        dis = sqrt(point(1,1)^2+point(1,2)^2+point(1,3)^2);
        if j==i+1
            %bounded calculations
            valVB = 0.5*k*(dis-le)^2;
            bGrad = -k*(dis-le);
            tempForceX = (bGrad*point(1,1)/dis);
            
            fX(i) = fX(i) + tempForceX;
            fX(j) = fX(j) - tempForceX;
        
            tempForceY = (bGrad*point(1,2)/dis);
            fY(i) = fY(i) + tempForceY;
            fY(j) = fY(j) - tempForceY;

            tempForceZ = (bGrad*point(1,3)/dis);
            fZ(i) = fZ(i) + tempForceZ;
            fZ(j) = fZ(j) - tempForceZ;
            tempP = tempP + valVB;
            
        else
            if ismember(i,pmatrix) || ismember(j,pmatrix)
                valV = 4*epsilonnb*((sigma/dis)^12-(sigma/dis)^6);
                val = 4*epsilonnb*((-12*sigma^12)/(dis^14)+(6*sigma^6)/(dis^8));
            else
                
                valV = 4*epsilon*((sigma/dis)^12-(sigma/dis)^6);
                val = 4*epsilon*((-12*sigma^12)/(dis^14)+(6*sigma^6)/(dis^8));
            end
            
            tempForceX = -(point(1,1)*val);
            fX(i) = fX(i) + tempForceX;
            fX(j) = fX(j) - tempForceX;

            tempForceY = -(point(1,2)*val);
            fY(i) = fY(i) + tempForceY;
            fY(j) = fY(j) - tempForceY;

            tempForceZ = -(point(1,3)*val);
            fZ(i) = fZ(i) + tempForceZ;
            fZ(j) = fZ(j) - tempForceZ;
        
            tempP = tempP + valV;
        end

    end
end