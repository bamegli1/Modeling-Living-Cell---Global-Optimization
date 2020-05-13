function[T] = trans_mat(states,stepsize,Nstructures,Kclusters)

T = zeros(Kclusters);

for k = 1:(Nstructures-stepsize)
    i = states(k);
    j = states(k+stepsize);
    T(i,j) = T(i,j) + 1;
end

T = T./sum(T(:));