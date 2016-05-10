function out = coagM0gainFrom0Alt(bins,alpha,d1,d2,rho,T,visc)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    
    if  d2/d1 > 3
        out = coagM0gainFrom0(bins,alpha,d1,d2,rho,T,visc);
        return;
    end
    
   
    if d2/d1-1<1e-3
        out=pi*(2*d1^3)*coag_kernel(pi/6*rho*d1^3,pi/6*rho*d2^3,d1,d2,T,visc);
        return;
    end
    
    c=alpha*ln(d2/d1);
    w=GaussOlinWeights4((1:4)',c);
    
    dVec=logspace(lg(d1),lg(d2),4);
    
    mp=pi/6*rho*dVec.^3;
    kern=coag_kernel(mp,mp,dVec,dVec,T,visc);
    
    [dp2,dp1]=meshgrid(dVec);
    kerroin=dp1.^3+dp2.^3;
    kern=kern.*kerroin;
    
    summa = w'*kern'*w;
    

%     for i=1:4
%         dp1=dVec(i);
%         
%         mp1=pi/6*rho*dp1^3;
%         
%         summa2=0;
%         for j=1:4
%             dp2=dVec(j);
%             
%             mp2=pi/6*rho*dp2^3;
%            
%             kern=coag_kernel(mp1,mp2,dp1,dp2,T,visc);
%             
%             summa2=summa2+(dp1^3+dp2^3)*kern*GaussOlinWeights4(j,c);
%             
%         end
%         
%         summa=summa+summa2*GaussOlinWeights4(i,c);
% 
%     end

    
    if abs(c)<1e-3
        out=pi*summa;
    else
        out = pi*summa*c^2/(exp(c)-1)^2;
    end

end
