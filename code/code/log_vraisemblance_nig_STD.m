%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
function [lv, g] = log_vraisemblance_nig_STD(x,data)
gama=exp(x(1));
betha=x(2);
delta=exp(x(3));
mu=x(4);
alpha=sqrt(gama^2+betha^2);
t = data-mu;
t=t';
t1=t/delta;
t2=t1.^2;
t3=sqrt(t2+1);
bes2=min(700,alpha*delta*t3);
bes2=max(exp(-30),bes2);

% Calcul des fonctions de Bessel
%-------------------------------
k2=besselk(1,bes2);
r2=besselk(2,bes2)./k2;
cst=pi;
n = length(data);
lv = -log(cst)+log(alpha)+delta*gama+mean(betha*t+log(k2)-log(t3));
if alpha+delta>50;
lv=-100;
end
%%%%/
% log(gama)=g_log
%-------
g(1) = gama*(2*gama/(alpha^2)+delta-(gama*delta/alpha)*mean(t3.*r2));
%%%%%
% betha
%-----
%
g(2)=2*(betha/(alpha^2))+mean(t-(betha*delta/alpha)*t3.*r2);
%%%%/
% log(delta)=d_log
%------
g(3)=delta*(gama+(1/delta)-(alpha)*mean(r2./t3));
%%%%%
% mu
%-------
g(4)=-betha+alpha*mean((t1./t3).*r2);

lv = -lv;
g = -g;

%disp(lv,'lv')
%disp(g,'g')
%disp([alpha betha delta mu],'[alpha betha delta mu]')

end
