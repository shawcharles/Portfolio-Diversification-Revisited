function [alpha,betha,delta,mu,m2,m1,m3,m4]=Normalisation_parametre(c,alpha0,betha0,delta0,mu0,info_betha)
% This transformation does not change the mathematical expectation, the variance, or the Skewness, but just the higher moments.

gama0=sqrt(alpha0^2-betha0^2);
%Calcul des trois 1er moments X de loi NIG(alpha0,betha0,delta0,mu0)
m2=(delta0*alpha0^2)/gama0^3;
m1=mu0+delta0*betha0/gama0;
m3=(3*delta0*betha0*alpha0^2)/gama0^5;
m4=3*delta0*alpha0^2*((alpha0^2+4*betha0^2)/gama0^7);

% Calculation of the parameters of the reduced centered variable Y=(X-m0)/s0
alpha=c*alpha0;

if info_betha==0
    if m3==0 
        betha=0;
    else
        Delta=(3*m2/m3)^2+4*alpha^2;
        b=3*m2/m3;
        betha_vect=.5*(-b+sqrt(Delta));
        betha_vect=[betha_vect,.5*(-b-sqrt(Delta))]
        diff=[alpha^2-betha_vect(1)^2,alpha^2-betha_vect(2)^2];
        [x,ind]=max(diff);
        betha=betha_vect(ind);
    end
else
    betha=c*info_betha;
end

gama=sqrt(alpha^2-betha^2);
delta=gama^3*(m2/alpha^2);
mu=m1-delta*betha/gama;

% We check that the first 3 moments have not varied.
m1=mu+delta*betha/gama
m2=(delta*alpha^2)/gama^3
m3=(3*delta*betha*alpha^2)/gama^5;
m4=3*delta*alpha^2*((alpha^2+4*betha^2)/gama^7);
m3=m3/m2^(3/2)
m4=m4/m2^2%excess Kurtosis

end
