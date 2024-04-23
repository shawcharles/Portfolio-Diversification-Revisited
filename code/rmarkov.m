function [x] = rmarkov(m,n)
% Simulates markov chain of length n with matrix m
x = ones(1,n);
for i = 2:n
   x(i) = rproba(m(x(i-1),:)); 
end
end