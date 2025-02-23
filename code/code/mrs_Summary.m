%=========================================================================
% Internally used routine(s)
%=========================================================================

%=========================================================================
% mrs_Summary
%=========================================================================

function Summary = mrs_Summary(Param,P,logL)

% Mean and variance
bParam = Param(1,:);
bPhi = bParam(1);
bC = bParam(2);
bSigma2 = bParam(3);
bMean =  bC / (1-bPhi);
bVar = nan;
sParam = Param(2,:);
sPhi = sParam(1);
sC = sParam(2);
sSigma2 = sParam(3);
sMean = sC / (1-sPhi);
sVar = nan;
% Unconditional probabilities
p_11 = P(1,1);
p_22 = P(2,2);
P_Rt_1 = (1-p_22) / (2-p_11-p_22);
P_Rt_2 = (1-p_11) / (2-p_11-p_22);
% Summary matrix
% Mean-reverting process
sPhi = sParam(1);
sC = sParam(2);
sSigma2 = sParam(3);        
Summary = [ bPhi,   bC,     bSigma2,  bMean,  bVar,   p_11,   P_Rt_1;
            sPhi,   sC,     sSigma2,  sMean,  sVar,   p_22,   P_Rt_2;
            logL,     nan,    nan,      nan,    nan,    nan,    nan];
end
