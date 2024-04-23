
function  [KS,Pvalue,tau, tau_rel] = TestChangePoint(X,N,est)
%
%  Nonparametric change point test for a univariate series  using the 
%  Kolmogorov-Smirnov  statistic.
%
% Input
%        X: (n x 1) vector of data (residuals or observations);
%        N: number of bootstrap samples to compute the P-value;
%      est: 1 if tau is estimated, 0 otherwise (default).
%
%  Output
%          KS: Kolmogorov-Smirnov statistic; 
%      Pvalue: (%) calculated with N bootstrap samples;
%         tau: estimation of change point time;
%     tau_rel: estimation of change point relative time (percent of the
%              sample size).
%
%   Ref: Remillard, B. (2012)  Non-Parametric Change Point Problems 
%        Using Multipliers, SSRN Working Paper Series No. 2043632.
%
%%

X=Data_final;
N=100
est=1
if nargin < 3
    est = 0;
end

n = length(X);


[KS,tau,tau_rel] = changepoint(X,est);

%%
KSsim = zeros(N,1);

    
    parfor k=1:N
        
        x = rand(n,1);
        tiedrank(x);
        KSsim(k) = changepoint(x,0);
    end
    
 Pvalue = 100*mean( (KSsim > KS) );
 
end
%%

function  [KS,tau,tau_rel] = changepoint(x,est)

n = length(x);
tau = [];
tau_rel = [];

R = tiedrank(x);

P = zeros(n,n);

z = P;

for k=1:n
    z(:,k) = (R <=k) -k/n;
end

P = cumsum(z); 
P = P/sqrt(n);  % process values!!

M1 = (max(abs(P')))';
KS= max(M1);

if(est)
tau=1;
while(M1(tau)< KS)
    tau = tau+1;
end

tau_rel = 100*tau/n;
end

end