function y = ig_rnd(n,delta,gama)
% y = ig_rnd(n,delta,gama) - 
% generation de v.a. inverse gaussienne
% la parametrisation de la densite est
% celle de la page 146 de la these de Raible
% ainsi que la methode de simulation
% attention pour la simulation il ya une erreure 
% dans la these de Raible (changer v-->v/2
% voir aussi p. 148 de Devroye

V=randn(1,n).^2;
if gama==0
  y=delta^2*ones(1,n)./V;
else
  U=rand(1,n);
  V=V/2;
  temp=V/(gama^2);
  z1=temp-sqrt(2*V*delta/gama^3+temp.^2)+delta/gama;
  z2=(delta/gama)^2*ones(1,n)./z1;
  p=delta*ones(1,n)./(gama*z1+delta);
  %disp([size(p),size(U)])
  I=find(p>U);
  y=z2;
  y(I)=z1(I);
end


