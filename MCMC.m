%Brandon Ameglio (bamegli1)
%Inputs = number of moves, kT value, box lengh, and a data set.
%Outputs = average potential energy and percentage acceptance of moves 
%MCMC function that generates moves and returns energies.

function[currV,currData]=MCMC(moves,kT,data)


pmatrix = [1,3,9,12,14,16,17,18,21,22,23,24,25,26];
sigma = 1;
epsilon = 1;

currData=data;

aveV = 0;
count = 0;

for n=1:moves
    %Random Particle and position
    randrow = randi(27);
    
    Delta = [(rand(1)-0.5)*0.1,(rand(1)-0.5)*0.1,...
        (rand(1)-0.5)*0.1];
    %Make movement
    newData = currData;
    newData(randrow,:) = newData(randrow,:)+Delta;
    
    %No PBC 
    %newData(randrow,:) = newData(randrow,:) - boxL * ...
    %    round(newData(randrow,:)/boxL);
    
    %Test energies    
    currV = energyFunctionMCMC(sigma,epsilon,currData,pmatrix);
    newV = energyFunctionMCMC(sigma,epsilon,newData,pmatrix);


    probACC = exp(-(newV-currV)/kT);
    if newV < currV | probACC >= (1)
        count = count + 1;
        currData = newData;
        aveV = aveV + newV;
    else
        aveV = aveV + currV;
    end
end
%count = (count/moves)*100;
%average = (aveV/moves);
