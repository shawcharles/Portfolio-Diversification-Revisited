function [a, lv, GRAD,std] = estim_nig(data)
% [alpha,betha,delta,mu], ecarts_types] = estim_nig(data)
%  estimer les parametres d'une distribution normale inverse gaussienne
% par la methode de maximum de vraisemblance
% la parametrisation de la densite est
% celle des lois hyper gene (voir these de Raible p. 144)
%
%%%%%%%%%%%%%%/
% Chgt de variable
% gama=sqrt(alpha^2-betha^2)
% on considere les para 
% gama, betha, delta, mu  
%%%%%%%%%%%%%%/
%INITIALISATION OPTIMISATION
%%%%%%%%%%%%%%/
%%%%%%%%%%%%%%%%%/
% methode=1 estimation des 4 para
% methode=2 betha est fixe=0
% methode=3 betha est fixe=0
%%%%%%%%%%%%%%%%%
methode=1;
%%%%%%%%%%%%%%%%/
%k2=(st_deviation(data))^2
%k1=mean(data);
%skew=mean((data-k1).^3)/(k2^(3/2))
%kurt=mean((data-k1).^4)/(k2^2)
%xsi=abs(3/(kurt-(4/3)*skew^2));
%ba=(skew^2)/(3*kurt-4*skew^2)
%%%%%%%%%%%%/
%alpha=3*sqrt(abs(ba))/abs((1-ba)*skew*sqrt(k2));
%gama=sqrt(abs(alpha^2-ba*alpha^2));
%betha=(skew*sqrt(k2)*gama^2)/3;
%delta=(gama^3/alpha^2)*k2;   
%mu=k1-delta*betha/gama    
%%%%%%%%%%%%%%/
mu1=mean(data);
mu2=mean((data-mu1).^2);
mu3=mean((data-mu1).^3);
mu4=mean((data-mu1).^4);
betha=1/((mu4/mu3)-5*mu3/(3*mu2));
gama=sqrt(abs(3*betha*mu2/mu3))
alpha=sqrt(betha^2+gama^2);
delta=(mu2*gama^3)/alpha^2;
mu=mu1-delta*betha/gama;

%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%E(data)=mu+delta*betha/gama;
%std(data)=sqrt((delta*alpha^2)/gama^3);
%skewness=3*delta*(alpha^2)*betha*gama^(-5)
%kurtosis=3*delta*(alpha^2)*(alpha^2+4*betha^2)*gama^(-7)
%skew=mean((data-k1).^3)
%kurt=mean((data-k1).^4)
%skweness2=3*alpha^(-1/4)*(betha/alpha)/(1-(betha/alpha)^2)^(1/4)
%kurtosis2=3*alpha^(-1/2)*(1+4*(betha/alpha)^2)/(1-(betha/alpha)^2)^(1/2)
%%%%%%%%%%%/
%delta1=sqrt(abs(k2*xsi*(1-rho^2)));
%alpha1=xsi/(delta*sqrt(abs(1-rho^2)));
%betha1=alpha1*rho;
%mu1=k1-rho*sqrt(abs(k2*xsi));
%disp([alpha betha delta mu],'alpha_init, betha_init, delta_init, mu_init')
%disp([alpha1 betha1 delta1 mu1],'alpha_init1, betha_init1, delta_init1, mu_init1')
%%%%%%%%%/
d_log=log(delta);
g_log=log(gama);

if methode==1
  %[lv ,aa,gradopt]=optim(list(log_vraisemblance_nig,data),[g_log betha d_log mu],'gc')
  [aa,lv,EXITFLAG,OUTPUT,GRAD,HESSIAN]=fminunc(@(x) Log_Vraisemb_Nig(x,data),[g_log betha d_log mu]);
  a(1)=sqrt(exp(aa(1))^2+aa(2)^2);
  a(2)=aa(2);
  a(3)=exp(aa(3));
  a(4)=aa(4);
elseif methode==2
  [lv ,aa,gradopt]=optim(list(log_vraisemblance_nig_betha,data),[g_log d_log mu],'gc')
  a(1)=exp(aa(1));
  a(2)=0;
  a(3)=exp(aa(2));
  a(4)=aa(3);
else
  %[lv ,aa,gradopt]=optim(list(log_vraisemblance_nig_betha2,data),[g_log d_log],'gc')
LB=[.01 .01];
UB=[20 20];
  [lv ,aa,gradopt]=optim(list(log_vraisemblance_nig_betha2,data),'b',LB,UB,[alpha delta],'gc')
  a(1)=aa(1);
  a(2)=0;
  a(3)=aa(2);
  a(4)=mean(data);
end
std=1;

end
