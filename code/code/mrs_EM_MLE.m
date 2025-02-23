%=========================================================================
% mrs_EM_MLE
%=========================================================================

function newParam = mrs_EM_MLE(Data,Ksi_tT,oldParam,gammas)
% MRS_EM_MLE ML estimation of regime processes parameters. 
sData = Data(1:end);
sKsi_tT = Ksi_tT(2:end, 2);

bData = Data(1:end);
bKsi_tT = Ksi_tT(2:end, 1);

if isempty(gammas),
    % Estimate g_1 and g_2
    % Spike regime
    sNewG = fminsearch(@(g) mle_g(sData,sKsi_tT,g),oldParam(2,4),optimset('Display','off'));
    sNewParam = mrEM_MLE_G(sData,sKsi_tT,sNewG);
    sNewPhi = 1-sNewParam(2);
    sNewC = sNewParam(1);
    sNewSigma2 = sNewParam(3);
    sNewParam = [sNewPhi, sNewC, sNewSigma2, sNewG,  nan];

    % Base regime
    bNewG = fminsearch(@(g) mle_g(bData,bKsi_tT,g),oldParam(1,4),optimset('Display','off'));
    bNewParam = mrEM_MLE_G(bData,bKsi_tT,bNewG);
    bNewPhi = 1-bNewParam(2);
    bNewC = bNewParam(1);
    bNewSigma2 = bNewParam(3);
    bNewParam = [bNewPhi, bNewC, bNewSigma2, bNewG, nan];
    
    %newParam = [bNewParam; sNewParam];
    newParam = [bNewPhi, bNewC, bNewSigma2, bNewG, nan;
        sNewPhi, sNewC, sNewSigma2, sNewG,  nan];
    
else
    % Use known g_1 and g_2
    % Spike regime
    sNewParam = mrEM_MLE_G(sData,sKsi_tT,oldParam(2,4));
    sNewPhi = 1-sNewParam(2);
    sNewC = sNewParam(1);
    sNewSigma2 = sNewParam(3);
    sNewParam = [sNewPhi, sNewC, sNewSigma2, oldParam(2,4), nan];

    % Base regime
    bNewParam = mrEM_MLE_G(bData,bKsi_tT,oldParam(1,4));
    bNewPhi = 1-bNewParam(2);
    bNewC = bNewParam(1);
    bNewSigma2 = bNewParam(3);
    bNewParam = [bNewPhi, bNewC, bNewSigma2, oldParam(1,4), nan];
    
    %newParam = [bNewParam; sNewParam];
    
    newParam = [bNewPhi, bNewC, bNewSigma2, oldParam(1,4), nan
        sNewPhi, sNewC, sNewSigma2, oldParam(2,4), nan];
    
end

end