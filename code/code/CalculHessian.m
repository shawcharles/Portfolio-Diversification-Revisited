function HESSIAN=CalculHessian(KSI_TT,Alpha,Betha,Sigma,lt,Delta,Data)
%Sigma est l'ecart type pas la variance
for K=1:2
    
    temp1=0;
    for t=2:lt
        temp1=temp1+KSI_TT(t,K)*Data(t-1)^(-2*Delta(K));
    end
    HESSIAN{K}(1,1)=-(1/(Sigma(K)^2))*temp1;
    
    temp2=0;
    for t=2:lt
        temp2=temp2+KSI_TT(t,K)*Data(t-1)^(1-2*Delta(K));
    end
    HESSIAN{K}(1,2)=(1/(Sigma(K)^2))*temp2;
    
    temp3=0;
    for t=2:lt
        temp3=temp3+KSI_TT(t,K)*(Data(t)-Alpha(K)-(1-Betha(K))*Data(t-1))*Data(t-1)^(-2*Delta(K));
    end
    HESSIAN{K}(1,3)=-(1/(Sigma(K)^4))*temp3;
        
    temp4=0;
    for t=2:lt
        temp4=temp4+KSI_TT(t,K)*Data(t-1)^(1-2*Delta(K));
    end
    HESSIAN{K}(2,1)=(1/(Sigma(K)^2))*temp4;
    
    temp5=0;
    for t=2:lt
        temp5=temp5+KSI_TT(t,K)*Data(t-1)^(2-2*Delta(K));
    end
    HESSIAN{K}(2,2)=-(1/(Sigma(K)^2))*temp5;
    
    temp6=0;
    for t=2:lt
        temp6=temp6+KSI_TT(t,K)*Data(t-1)^(-2*Delta(K))*((Data(t)-Alpha(K)-(1-Betha(K))*Data(t-1))*Data(t-1));
    end
    HESSIAN{K}(2,3)=(1/(Sigma(K)^4))*temp6;
    
    temp7=0;
    for t=2:lt
        temp7=temp7+KSI_TT(t,K)*Data(t-1)^(-2*Delta(K))*(Data(t)-Alpha(K)-(1-Betha(K))*Data(t-1));
    end
    HESSIAN{K}(3,1)=-(1/(Sigma(K)^4))*temp7;
    
    temp8=0;
    for t=2:lt
        temp8=temp8+KSI_TT(t,K)*Data(t-1)^(1-2*Delta(K))*(Data(t)-Alpha(K)-(1-Betha(K))*Data(t-1));
    end
    HESSIAN{K}(3,2)=(1/(Sigma(K)^4))*temp8;
    
    temp9=0;
    temp10=0;
    for t=2:lt
        temp10=temp10+KSI_TT(t,K);
        temp9=temp9+KSI_TT(t,K)*Data(t-1)^(-2*Delta(K))*(Data(t)-Alpha(K)-(1-Betha(K))*Data(t-1))^2;
    end
    HESSIAN{K}(3,3)=(1/(2*Sigma(K)^4))*temp10-(1/(Sigma(K)^6))*temp9;
end
end