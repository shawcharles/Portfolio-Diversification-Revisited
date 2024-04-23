clc
clear all
close all
cd 'D:\AVINASH_PAPER_CODE\Chapter 2\code'
graine=888;
rand('seed',graine);
randn('seed',graine);

MC=10;
ppp=0.1;


%X=xlsread('mordata.xls','BB','b797:b5810');
X=xlsread('newdata.csv','newdata','b2:b5164');
%X=100*(log(equity(2:end,:))-log(equity(1:end-1,:)));
%Data=(log(equity(2:end,:))-log(equity(1:end-1,:)));



pas=1;
Dt=1/252;
DataOut=CutData(X,pas);
Data=DataOut';

DeltaChoix=4;

% Estimation of the RS Dynamic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[KSI_TT,PARAM,P,KSI_T1T_10,LOGL,KSI_tt,PARAM_stock]=ps2r_est(Data,[0 0]);

%Smoothed probability: KSI_TT
figure(1)
subplot(2,1,1)

hold on
grid
title('Filtered Probability')
plot(KSI_tt(:,1),'LineWidth',2);
plot(KSI_tt(:,2),'r','LineWidth',2);
legend('Regime 1','Regime 2');
subplot(2,1,2)
hold on
grid

title('Smoothed Probability')
plot(KSI_TT(:,1),'LineWidth',2);
plot(KSI_TT(:,2),'r','LineWidth',2);
legend('Regime 1','Regime 2');

Betha(1)=1-PARAM(1,1);
Betha(2)=1-PARAM(2,1);
% Unconditional probabilities
p_11 = P(1,1);
p_22 = P(2,2);
P_Rt_1 = (1-p_22) / (2-p_11-p_22);
P_Rt_2 = (1-p_11) / (2-p_11-p_22);

Mean(1) =  PARAM(1,2) / (1-PARAM(1,1));
Mean(2) =  PARAM(2,2) / (1-PARAM(2,1));



%Calcul of the Std
DeltaModele{1}=[PARAM(1,4) PARAM(2,4)];
AlphaModele{1}=[PARAM(1,2) PARAM(2,2)];
BethaModele{1}=[Betha(1) Betha(2)];
SigmaModele{1}=[sqrt(PARAM(1,3)) sqrt(PARAM(2,3))];
lt=length(KSI_TT(:,1));

%Evaluation of the Hessian Matrix
HESSIAN=CalculHessian(KSI_TT,AlphaModele{1},BethaModele{1},SigmaModele{1},lt,DeltaModele{1},Data);

%Fisher Information Matrix
Fisher{1}=-HESSIAN{1};
Fisher{2}=-HESSIAN{2};
STD_M{1}=sqrt(inv(Fisher{1}));
STD_ERROR{1}{1}=[STD_M{1}(1,1) STD_M{1}(2,2) STD_M{1}(3,3)];
STD_M{2}=sqrt(inv(Fisher{2}));
STD_ERROR{1}{2}=[STD_M{2}(1,1) STD_M{2}(2,2) STD_M{2}(3,3)];
STD_MODELE=[STD_M{1}(1,1) STD_M{1}(2,2) STD_M{1}(3,3);STD_M{2}(1,1) STD_M{2}(2,2) STD_M{2}(3,3)];



fprintf('Parameters of our model\n')
fprintf('               Regime 1    Regime 2\n')
fprintf('Alpha:         %f          %f      \n',PARAM(1,2),PARAM(2,2));
fprintf('Betha:         %f          %f      \n',Betha(1),Betha(2));
fprintf('Sigma:         %f          %f      \n',sqrt(PARAM(1,3)),sqrt(PARAM(2,3)));
fprintf('Et             %f          %f      \n',Mean(1),Mean(2));
fprintf('q_ii           %f          %f      \n',p_11,p_22);
fprintf('P(R=i)         %f          %f      \n',P_Rt_1,P_Rt_2);


%Evolution of Parameters
NN=length(PARAM_stock);
for p=1:NN
    AlphaTheta(:,p)=PARAM_stock{p}(:,2);
    BethaTheta(:,p)=1-PARAM_stock{p}(:,1);
    SigmaTheta(:,p)=sqrt(PARAM_stock{p}(:,3));
end
figure(2)
subplot(3,1,1)
hold on
grid
title('Alpha')
plot(AlphaTheta(1,:),'LineWidth',2);
plot(AlphaTheta(2,:),'r','LineWidth',2);
legend('Regime 1','Regime 2');
subplot(3,1,2)
hold on
grid
title('Beta')
plot(BethaTheta(1,:),'LineWidth',2);
plot(BethaTheta(2,:),'r','LineWidth',2);
legend('Regime 1','Regime 2');
subplot(3,1,3)
hold on
grid
title('Sigma')
plot(SigmaTheta(1,:),'LineWidth',2);
plot(SigmaTheta(2,:),'r','LineWidth',2);
legend('Regime 1','Regime 2');
m1=length(Data);
t1=1:1:m1;
Date=1:36:m1;
Date=[Date m1];
figure(3)
hold on;
grid
plot(t1,Data,'LineWidth',2);
xlabel('Dates')


%Graph of the State
%%%%%%%%%%%%%%%%%%%%
ll=length(KSI_TT(:,1));
for q=1:ll
    Ind(q)=find(KSI_TT(q,:)==max(KSI_TT(q,:)));
end
figure(4)
hold on
grid
plot(Ind,'+');

figure(5)
hold on
grid
DataState1=[];
DataState2=[];
for k=q:ll
    if Ind(q)==1
        plot(q,Data(q),'+','LineWidth',2);
        DataState1=[DataState1;Data(q)];
    else
        plot(q,Data(q),'r+','LineWidth',2);
        DataState2=[DataState2;Data(q)];
    end
end

xlabel('Dates')



figure(6)
subplot(3,1,1)
hold on
grid
title('Filtered Probability')
plot(KSI_tt(:,1),'LineWidth',2);
plot(KSI_tt(:,2),'r','LineWidth',2);
legend('Regime 1','Regime 2');
subplot(3,1,2)
hold on
grid
title('Smoothed Probability')
plot(KSI_TT(:,1),'LineWidth',2);
plot(KSI_TT(:,2),'r','LineWidth',2);
legend('Regime 1','Regime 2');
subplot(3,1,3)
hold on
grid
DataState1=[];
DataState2=[];
for q=1:ll
    if Ind(q)==1
        plot(q,Data(q),'+','LineWidth',2);
        DataState1=[DataState1;Data(q)];
    else
        plot(q,Data(q),'r+','LineWidth',2);
        DataState2=[DataState2;Data(q)];
    end
end
xlabel('Dates')




%Long Run mean
fprintf('Theta:         %f         %f\n',PARAM(1,2)/Betha(1),PARAM(2,2)/Betha(2));
fprintf('Kappa:         %f         %f\n',Betha(1),Betha(2));
fprintf('Sigma:         %f         %f\n',PARAM(1,3),PARAM(2,3));

RESPARAM=[Betha(1) PARAM(1,2)/Betha(1) PARAM(1,3) p_11; Betha(2) PARAM(2,2)/Betha(2) PARAM(2,3) p_22];


%Calcul RCM statistic
%%%%%%%%%%%%%%%%%%%%%
lt=length(KSI_TT(:,1));
K=2;
RCM=0;
for q=1:lt
    RCM=RCM+KSI_TT(q,1)*KSI_TT(q,2);
end
RCM=RCM*(100*2^2)/(lt);
RCM2bis=0;
for q=1:lt
    RCM2bis=RCM2bis+(KSI_TT(q,1)-1/K)^2+(KSI_TT(q,2)-1/K)^2;
end
RCM2=100*(1-((K)/(K-1))*(1/lt)*RCM2bis);
%fprintf('RCM Statistic 1: %f \n',RCM);
fprintf('RCM Statistic: %f \n',RCM2);

Perct=0;

for q=1:lt
    for l=1:K
        if KSI_TT(q,l) > 1-ppp
            Perct=Perct+1;
        elseif KSI_TT(q,l) <ppp
            Perct=Perct+1;
        end
    end
end
PerctFinal=(Perct/(lt*K))*100;
fprintf('Percentage of classification: %f\n',PerctFinal);


%TABLE 1 STATISTIC GENERAL
x1=Data;
[nobs1,junk] = size(x1);
maxim=  max(x1);
minim = min(x1);
mean1 = 1/nobs1*sum(x1);
var1 = 1/nobs1*sum((x1-mean1).^2);
stddev1 = sqrt(var1);
skew = mean((x1-mean1).^3);
kurtosis = mean((x1-mean1).^4);

fprintf('Total number of observations: %f\n',nobs1);
fprintf('maximum: %f\n',maxim);
fprintf('minimum: %f\n',minim);
fprintf('mean: %f\n',mean1);
fprintf('variance: %f\n',var1);
fprintf('standard deviation: %f\n',stddev1);
fprintf('Skewness: %f\n',skew);
fprintf('Kurtosis: %f\n',kurtosis);

%TestLevy=input('\n Test de la présence de sauts ? \n 1- Oui \n 2 - Non \n'),

%Levy Estimation
%%%%%%%%%%%%%%%%%%%%%%%%
TestLevy=1;
if TestLevy==1
    
    for ll=1:2;
               %Step 1: Least squared method
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ll==1
            X=DataState1;
        else
            X=DataState2;
        end
        Length=length(X);
        A=[ones(Length-1,1) X(end-1:-1:1)];
        B=X(end:-1:2);
        Res=inv(transpose(A)*A)*transpose(A)*B;
        hat_m=Res(1);
        hat_a=1-Res(2);
        Temp=(X(2:end)-(1-hat_a).*X(1:end-1)).^2;
        Temp2=cumsum(Temp);
        hat_s2=Temp2(end)/Length;
        %Step 2: MLE
        %%%%%%%%%%%%
        %Construction of length Y sequence 
        Y=X(2:end)-(1-hat_a)*X(1:end-1);
        %Paramètres Initiaux
        Y0=cumsum(Y);
        bar_Y=Y0(end)/Length;
        Y1=cumsum(Y-bar_Y);
        Y2=cumsum((Y-bar_Y).^2);
        Y3=cumsum((Y-bar_Y).^3);
        Y4=cumsum((Y-bar_Y).^4);
        YY1=Y1(end)/Length;
        YY2=Y2(end)/Length;
        YY3=Y3(end)/Length;
        YY4=Y4(end)/Length;
        bar_S=YY2;
        bar_YYY1=YY3/(YY2^(3/2));
        bar_YYY2=(YY4/(YY2^2))-2;
        
        gama=3/(bar_S*sqrt(3*bar_YYY2-5*bar_YYY1^2));
        betha=(bar_YYY1*bar_S*gama^2)/(3);
        delta=(bar_S^2*gama^3)/(betha^2+gama^2);
        mu=bar_Y-betha*((delta)/(gama));
        alpha=sqrt(betha^2+gama^2);
        d_log=log(delta);
        g_log=log(gama);
        lam=1;
        
        %Estimation
        param=[g_log betha d_log mu];
        options = optimoptions(@fminunc,'Display','off');
        %[est,FVAL,EXITFLAG,OUTPUT,GRAD,HESSIAN]=fminunc(@(x) Log_Vraisemb_Nig(x,Y),param);
        [est,FVAL]=fminunc(@(x) Log_Vraisemb_Nig(x,Y),param,options);
        
        a(ll,1)=sqrt(exp(est(1))^2+est(2)^2);
        a(ll,2)=est(2);
        a(ll,3)=exp(est(3));
        a(ll,4)=est(4);
        hat_s=sqrt(hat_s2);
        NIG(ll,1)=hat_s*a(ll,1);
        NIG(ll,2)=hat_s*a(ll,2);
        NIG(ll,3)=a(ll,3)/hat_s;
        NIG(ll,4)=(a(ll,4)-hat_m)/hat_s;
        NIG(ll,5)=sqrt(NIG(ll,1)^2-NIG(ll,2)^2);
        
        
        [pdfNIG,g]=Log_Vraisemb_Nig([log(sqrt(NIG(ll,1)^2-NIG(ll,2)^2)),NIG(ll,2),log(NIG(ll,3)),NIG(ll,4)],Y);
        LLog(ll)=pdfNIG;
        
        kappa(ll)=-(1/Dt)*log(1-hat_a);
        theta(ll)=hat_m/hat_a;
        sigma(ll)=hat_s*(1/sqrt(hat_a/(2*kappa(ll))));
        
        
        gama(ll)=sqrt(NIG(ll,1)^2+NIG(ll,2)^2);
        d_log(ll)=log(NIG(ll,3));
        g_log(ll)=log(gama(ll));
		        for compt=1:length(Y);
            [lv, g] = log_vraisemblance_nig_STD([g_log(ll) NIG(ll,2) d_log(ll) NIG(ll,4)],Y(compt));
            Hess_STD(compt,:)=g;
        end
        Hess_Mat{ll}=zeros(4,4);
        for compt=1:length(Y);
            Hess_Mat{ll}=Hess_Mat{ll}+Hess_STD(compt,:)'*Hess_STD(compt,:);
        end
        Hess_Mat{ll}=Hess_Mat{ll}./length(Y);                
    end
    
		
    fprintf('Parameters of our model NIG\n')
    fprintf('               Regime 1    Regime 2\n')
    fprintf('Alpha             %f          %f      \n',NIG(1,1),NIG(2,1));
    fprintf('Betha             %f          %f      \n',NIG(1,2),NIG(2,2));
    fprintf('Delta             %f          %f      \n',NIG(1,3),NIG(2,3));
    fprintf('Mu                %f          %f      \n',NIG(1,4),NIG(2,4));
    
    
    %TABLE 18 AIC BIC Results
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    LOG_LEVY=-(exp(LLog(1))+exp(LLog(2)));
    fprintf('\nLog Likelihood RS-Levy: %f\n',LOG_LEVY);
    AIC_RS_L=-2*LOG_LEVY+2*8;
    BIC_RS_L=-2*LOG_LEVY+8*log(length(Data));
    
    AIC_RS_G=-2*LOGL+2*8;
    BIC_RS_G=-2*LOGL+8*log(length(Data));
    
   
    fprintf('RS_Levy Case BIC: %f and AIC %f \n',BIC_RS_L,AIC_RS_L);
    fprintf('RS_Gaussian Case BIC: %f and AIC %f \n',BIC_RS_G,AIC_RS_G);

    %Gaussian Case
    DELTAT=1;
    [mu2,sigma2,lambda2]=OU_Calibrate_ML(Data',DELTAT);
    sigma2_bar=sigma2^2*((1-exp(-2*lambda2*DELTAT))/(2*lambda2));
    SUM_V=0;
    for i=2:lt;
        SUM_V=SUM_V+(Data(i)-Data(i-1)*exp(-lambda2*DELTAT)-mu2*(1-exp(-lambda2*DELTAT)))^2;
    end
    logV_vasicek=-lt*log(sqrt(sigma2_bar))-(lt/2)*log(2*pi)-(1/(2*sigma2_bar))*SUM_V;
    a=lambda2*mu2;
    b=lambda2;
    fprintf('Log Likelihood Gaussian Case: %f\n',logV_vasicek);    
    AIC_G=-2*logV_vasicek+2*3;
    BIC_G=-2*logV_vasicek+3*log(length(Data));
    fprintf('Gaussian Case BIC: %f and AIC %f \n',BIC_G,AIC_G);         
end

Hess_Mat_Final=[Hess_Mat{1}(1,1) Hess_Mat{2}(1,1);Hess_Mat{1}(2,2) Hess_Mat{2}(2,2);Hess_Mat{1}(3,3) Hess_Mat{2}(3,3);Hess_Mat{1}(4,4) Hess_Mat{2}(4,4)];

diary myDiaryFile