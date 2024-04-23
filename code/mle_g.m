%==========================================================================
% mrs_EM_MLE auxiliary functions 
%==========================================================================

function LogL = mle_g(Data,KsiT_t,g)
beta = (sum(KsiT_t.*Data(1:end-1).*abs(Data(1:end-1)).^(-2*g).*( Data(2:end) - Data(1:end-1))) -...
    sum( KsiT_t ./ abs(Data(1:end-1)).^(2*g) .* ( Data(2:end) - Data(1:end-1) ) ) ./ sum ( KsiT_t ./ abs(Data(1:end-1)).^(2*g) ).*...
    sum(KsiT_t .*Data(1:end-1).*abs(Data(1:end-1)).^(-2*g) ))./...
    (sum(KsiT_t.*Data(1:end-1).*abs(Data(1:end-1)).^(-2*g)).* sum(KsiT_t.*Data(1:end-1).*abs(Data(1:end-1)).^(-2*g)) ./ sum( KsiT_t ./ abs(Data(1:end-1)).^(2*g) ) - sum(abs(Data(1:end-1)).^(2-2*g).*KsiT_t)  );
alpha = sum( KsiT_t ./ abs(Data(1:end-1)).^(2*g) .* ( Data(2:end) - Data(1:end-1) ))./sum(KsiT_t./abs(Data(1:end-1)).^(2*g))...
    + beta * sum(Data(1:end-1).*abs(Data(1:end-1)).^(-2*g).*KsiT_t )  ./ sum(KsiT_t./abs(Data(1:end-1)).^(2*g));
sigma2 = sum(KsiT_t./abs(Data(1:end-1)).^(2*g).*(Data(2:end)-Data(1:end-1)- alpha + beta*Data(1:end-1)).^2)./sum(KsiT_t);
LogL = abs(-sum(KsiT_t.*log(abs(Data(1:end-1))))+sum(log(abs(Data(1:end-1))).*(Data(2:end)-(1-beta).*Data(1:end-1)-alpha).^2./sigma2./abs(Data(1:end-1)).^(2*g).*KsiT_t));
end
