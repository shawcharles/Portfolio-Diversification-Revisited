function x_simu = simu_OU_NIG_RS(n,sigma,mu,Z)
% simulation of a Brownian process
% X(t+1)=X(t)+mu+sigma*eps(t)
% eps_t gaussian centered reduced iid
% n=growth simulator
% sigma = volatility (possibly multidimensional, if periodic)
% mu = average (possibly multidimensional, if periodic)
eps=randn(1,n-1);
eps=nigrnd(alpha,beta,mu,delta,1,n-1)

for k=1:n-1
    if Z==1
        mu_RS(k)=mu(1);
        sigma_RS(k)=sigma(1);
    else
        mu_RS(k)=mu(2);
        sigma_RS(k)=sigma(2);
    end
end
x_simu=mu_RS.*(1:n-1)+sigma_RS.*cumsum(eps);

x_simu=[0, x_simu];
end

