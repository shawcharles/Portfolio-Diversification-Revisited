function [x] = Simul_markov(m,n,i)
% Simulates a markov chain of length n with matrix m and initial state i
x = ones(1,n);
x(1)=i;
for i = 2:n
   x(i) = rproba(m(x(i-1),:)); 
end
end