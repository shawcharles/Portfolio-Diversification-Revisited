function [llrstat, llrprb] = llrtest(llr,kr,llu,ku)
% PURPOSE: computes likelihood ratio test for two regressions
%---------------------------------------------------
% USAGE: [llrstat llrprob] = llr(llr,kr,llu,ku)
%    or: llr(llr,kr,llu,ku), which prints to the screen
% Where: llr = negative of the log likelihood from restricted model
%        kr = number of parameters in the restricted model
%        llu = negative of the log likelihood from unrestricted model
%        ku = number of parameters in the unrestricted model
%---------------------------------------------------
% RETURNS: llrstat = {-2(-logLrestricted+logLunrestricted)} 
%          llrprb  = marginal probability for llrstat
% NOTE:  large llrstat => reject the restrictions as inconsisent 
%                       with the data
% NOTE:  usually this is presented as
%        -2(logLrestricted-logLunrestricted), but llrtest requires
%        the negative of the log likelihood
% NOTE: requires chis_prb from Jim Lesage's Econometrics toolbox
%---------------------------------------------------

% written by:
% Guillaume Frechette
% Ohio State University
% Department of Economics
% 410 Arps Hall
% 1945 North High Street
% Columbus, OH 43210-1172
% frechette.6@osu.edu

pflag = 0;
if nargout == 0
pflag = 1;
else nargin ~= 4 % flag incorrect arguments
    error('llrtest: Wrong # of input arguments');
end;
if (ku - kr) < 0 % flag reversed input arguments
error('llrtest: negative dof, check for reversed input arguments');
end;
% check that llr>=llu (don't forget ll really is the negative of the loglikelihood
if (llr-llu)<0, error('llrtets: restricted likelihood is greater than the unrestricted one');end;
% compute the llrstat
numr = ku - kr; % find # of restrictions
llrstat = -2 * ( -llr + llu ); llrprb = chis_prb(llrstat,numr);

if pflag == 1
fprintf(1,'log likelihood ratio stat = %16.8f \n',llrstat);
fprintf(1,'probability               =%16.4f\n',llrprb);
fprintf(1,'num, dof                  =         %4d\n',numr);
end;
