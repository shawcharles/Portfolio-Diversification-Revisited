function F=CalculPDFCIR(Theta,Datat1,Datat)

%Nombre de regime
[m,n]=size(Theta);

alpha=Theta(:,1);
betha=Theta(:,2);
sigma=Theta(:,3);


for i=1:m
   AA=sqrt(2*pi)*sigma(i)*abs(Datat)^.5;
   B1=Datat1-(1-betha(i))*Datat-alpha(i);
   B2=2*sigma(i)^2*abs(Datat);
   F(i)=(1/AA)*exp(-(B1)^2/(B2));
end

end
