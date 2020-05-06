%Brandon Ameglio (bamegli1)
%Inputs = sigma, epsilon, box lengh, and a data set.
%Outputs = potential energy of data. 
%Energy funciton for finding LJP of a dataset based on distances. 
%LD energy function of a particle

function[tempP]=energyFunctionMCMC(sigma,epsilon,data,pmatrix)

k = 20;
le = 1;
epsilonnb = 2/3;

tempP = 0;
for i=1:length(data)
    for j=i+1:length(data)
        point = [data(i,1)-data(j,1),data(i,2)-data(j,2),data(i,3)-data(j,3)];
        dis = sqrt(point(1,1)^2+point(1,2)^2+point(1,3)^2);
        if j==i+1
            %bounded calculations
            valVB = 0.5*k*(dis-le)^2;          
            tempP = tempP + valVB;
            
        else
            if ismember(i,pmatrix) || ismember(j,pmatrix)
                valV = 4*epsilonnb*((sigma/dis)^12-(sigma/dis)^6);
            else
                valV = 4*epsilon*((sigma/dis)^12-(sigma/dis)^6);
            end      
            tempP = tempP + valV;
        end

    end
end