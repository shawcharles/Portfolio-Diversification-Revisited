function [Ksi_tT,Param,P,Ksi_t1t_10,LogL,Ksi_tt,Param_stock]=ps2r_est(Data,gammas,display,startParam,startP,startKsi_t1t_10)
%PS2R_EST Estimates parameters of a 2-regime parameter switching (PS) model. 
%   [KSI_TT,PARAM,P]=PS2R_EST(DATA) returns smoothed inferences KSI_TT, 
%   estimated parameters PARAM and transition matrix P of a 2-regime 
%   parameter switching (PS) model, i.e. a 2-regime Markov regime-switching 
%   (MRS) model with both regimes driven by AR(1) processes of the form: 
%   X(t+1)=phi_i*X(t)+c_i+sigma_i*|X(t)|^g_i*N(0,1).
%   The first column (KSI_TT) or row (PARAM, P) contains results for the
%   base regime and the second column/row for the spike regime.
%   [KSI_TT,PARAM,P,KSI_T1T_10,LOGL]=PS2R_EST(DATA) additionally returns 
%   probabilities KSI_T1T_10 classifying the first observation to one of 
%   the regimes and log-likelihood LOGL of the fitted model. 
%   
%   [...]=PS2R_EST(DATA,GAMMAS) allows to specify the known values of g_1
%   and g_2 (as 1x2 vector [g_1, g_2]; default value GAMMAS=[]). If vector 
%   GAMMAS is not given (or GAMMAS=[]) g_1 and g_2 are estimated, otherwise 
%   a numerically less demanding algorithm is used for estimating the 
%   remaining parameters.
%
%   PS2R_EST(DATA,GAMMAS,DISPLAY,STARTPARAM,STARTP,STARTKSIT1T_10) allows 
%   to specify initial (model-dependent) parameter estimates STARTPARAM 
%   (a 2x5 vector), initial transition matrix STARTP and initial estimates
%   STARTKSIT1T_10 of probabilities classifying the first observation.
%   Default values for the latter three parameters are:
%       STARTPARAM = [0.3, 15, 1, 0, NaN; 0.3, 15, 1, 0, NaN];
%       STARTP = [0.8, 0.2; 0.6, 0.4];
%       STARTKSIT1T_10 = [0.6, 0.4];
%   The third pameter, DISPLAY (default DISPLAY=1), is a flag which defines 
%   whether the calibration results are dispayed in the command window (1) 
%   or not (0).
%   
%   Example:
%       Param = [0.2,2,1,0;0.4,3,1,1]; P=[0.5,0.5;0.4,0.6];
%       [Y,S] = ps2r_sim(P,Param,3.5,[1,0],1000,1);
%       ps2r_est(Y);
%       ps2r_est(Y,[0,1]);  
%   
%   See also PS2R_SIM
%
%   Reference(s):
%   [1] J.Janczura, R.Weron (2011) Efficient estimation of Markov 
%   regime-switching models: An application to electricity spot prices. 
%   Working paper version available at: 
%   http://ideas.repec.org/p/wuu/wpaper/hsc1102.html

%   Written by Joanna Janczura and Rafal Weron (2011.02.28)
%   Revised by Joanna Janczura and Rafal Weron (2011.10.03)

if nargin < 6
    startKsi_t1t_10 = [0.6, 0.4];
end
if nargin < 5 
    startP = [0.5, 0.5; 0.2, 0.8];
end
if nargin<2
    gammas = [];
end
if nargin < 4
    startParam = [0.3, 15, 1, 0, NaN; 0.3, 15, 1, 0, NaN];
    T = mean(Data);
    s = std(diff(Data));
    Pind = (Data>(T+s));
    startParam = mrs_EM_MLE(Data, [1-Pind,Pind], startParam, gammas);
end
if nargin<3
    display = 1;
end
if numel(gammas)==2,
    % Assign known values of g_1 and g_2
    startParam(1,4) = gammas(1);
    startParam(2,4) = gammas(2);
else
    % Estimate g_1 and g_2
    gammas = [];
end

% Initial estimation step 
[newKsi_tT,  newP, newKsi_t1t_10, newLogL,newKsi_tt] = mrs_EM_Smoother(Data, startParam, startP, startKsi_t1t_10);   
newParam = mrs_EM_MLE(Data, newKsi_tT, startParam, gammas);

Param_stock{1}=newParam;

% Estimation loop 
iteration = 1;
dif = 1;
KK=100;
while and(dif > 10^(-8), iteration < KK)
    oldP = newP;
    oldKsi_t1t_10 = newKsi_t1t_10;
    oldParam = newParam;

    oldKsi_tT = newKsi_tT;
    oldLogL = newLogL; 
    [newKsi_tT, newP, newKsi_t1t_10, newLogL,newKsi_tt] = mrs_EM_Smoother(Data, oldParam, oldP, oldKsi_t1t_10);    
    newParam = mrs_EM_MLE(Data, newKsi_tT, oldParam, gammas);
    
    %newParam(1,2)=newParam(1,1);
    
    
    %Stockage des paramètres
    Param_stock{iteration+1}=newParam;
    
    diffP = max(abs(newP-oldP));
    diffKsi_t1t_10 = max(abs(newKsi_t1t_10-oldKsi_t1t_10));
    diffParam = max(abs(newParam-oldParam));
    diffKsi_tT = max(abs(newKsi_tT-oldKsi_tT));
    diffLogL = abs(newLogL-oldLogL);
    dif = max([diffP diffKsi_t1t_10 diffParam diffKsi_tT diffLogL]);
    iteration = iteration+1;
    

    
end

% Results
Ksi_tt = newKsi_tt;
Ksi_tT = newKsi_tT;
Param = newParam;
P = newP;
ind = find(P<0);
if ~isempty(ind)
    P(ind) = 0;
end
Ksi_t1t_10 = newKsi_t1t_10;
LogL = newLogL;
% Display results in the command window
if display
    Summary = mrs_Summary(Param, P, LogL);
    disp([' ']);
    disp(['Gamma=' sprintf('%0.2g',Param(1,4)) ' for base regime']);
    disp(['Gamma=' sprintf('%0.2g',Param(2,4)) ' for spike regime']);
    % mean-reverting process
    disp([' ']);        
    disp(['Two state regime switching model with mean-reverting process for spikes, LogL=' num2str(LogL)]);                        
    disp([' ']);        
    disp(['regime   phi_i   c_i      sigma^2_i  E(Y_{t,i})  Var(Y_{t,i})  q_ii     P(R=i)']);
    disp(['base    ' sprintf('%7.5f',Summary(1,1)) '  ' sprintf('%7.5f',Summary(1,2)) '  ' sprintf('%7.5f',Summary(1,3)) '    ' ...
        sprintf('%7.5f',Summary(1,4)) '     ' sprintf('%7.5f',Summary(1,5)) '       ' sprintf('%7.5f',Summary(1,6)) '  ' sprintf('%7.5f',Summary(1,7))]);
    disp(['spike   ' sprintf('%7.5f',Summary(2,1)) '  ' sprintf('%7.5f',Summary(2,2)) '  ' sprintf('%7.5f',Summary(2,3)) '    ' ...
        sprintf('%7.5f',Summary(2,4)) '     ' sprintf('%7.5f',Summary(2,5)) '       ' sprintf('%7.5f',Summary(2,6)) '  ' sprintf('%7.5f',Summary(2,7))]);    
end
end

