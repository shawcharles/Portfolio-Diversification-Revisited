function [x] = rproba(p)
% draws a sample according to the law p carried by [1, N]
u = rand();
cdf = cumsum(p);
x = sum(cdf<u)+1;
end