function x = nig_rnd(n, alpha,bbeta,mu,delta)
% x = nig_rnd(n,alpha,beta,mu,delta)
% Normal Gaussian inverse with parametrisation
% Generate a sequence of n v.a distributed according to the law
% variance gamma N(mu+beta z,z) or z follows a Gaussian inverse


gama=sqrt(alpha^2-bbeta^2);
z=ig_rnd(n,delta,gama);
%disp(sum(imag(z)),'sum(imag(z))')
x=bbeta*z+mu+sqrt(z).*randn(1,n);

%/////////////////////////
% test
%/////////////////////////
%stacksize(30000000);getf("/home/oudjane/SCILAB/MODEL/SPOT/nig_rnd.sci"); getf("/home/oudjane/SCILAB/MODEL/SPOT/nig_pdf.sci"); getf("/home/oudjane/SCILAB/MODEL/SPOT/kernel.sci"); n=1000;  alpha=.5; beta=0;  mu=1; delta=1; x = nig_rnd(n, alpha,beta,mu,delta); pas=max(x)/1000; y=min(x):pas:max(x); z=kernel(x',y');zz=nig_pdf(y,alpha,beta,mu,delta); xset("window",2); plot2d(y,zz);plot2d(y,z,-2);  





