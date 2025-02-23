y=xlsread('E:\RSPAPER\CODE\data_c.csv')

% do test
[h,pValue] = adftest(y)

% compute log returns, test again then add lags
r = diff(log(y))

[h,pValue] = adftest(r)

% Graph for Markov Switching model
y=xlsread('E:\RSPAPER\CODE\data_c.csv')
returns = diff(log(y))

indep = ones(size(returns)); % A dummy explanatory variable
k = 2; % How many regimes we expect: bull and bear
S = [1 1]; % Both the mean and the volatility differ in bulls and bears

SpecOut = MS_Regress_Fit(returns, indep, k, S);