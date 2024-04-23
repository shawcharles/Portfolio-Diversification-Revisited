% Helper functions

function output = logdensity2info(logdensity,x,del,param)
% Compute the information matrix = Sum_{i=1}^n Score_i' Score_i
n = length(x) - 1;
output = 0;
for i=1:n
    tmpfun = @(theta)(logdensity(x(i+1,:),x(i,:),del,theta));
    score_i = gradest(tmpfun,param); % This is a row vector
    output = output + score_i' * score_i;
end