function [variance] = variancefun(V)


n = length(V);
m = mean(V);
x = 0;

for i = 1:n
    x = x +(V(i)-m)^2;
end

variance = x/n;

