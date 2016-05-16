function out = coagN0Alt(bins,alpha,d1,d2,rho,T,visc)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    if  d2/d1 > 3
        out = coagN0(bins,alpha,d1,d2,rho,T,visc);
        return;
    end
   
    if d2/d1-1<1e-3
        out=-pi*coag_kernel(pi/6*rho*d1^3,pi/6*rho*d2^3,d1,d2,T,visc);
        return;
    end
    
    c=alpha*ln(d2/d1);
    w=GaussOlinWeights4((1:4)',c);
    
    dVec=logspace(lg(d1),lg(d2),4);
    
    mp=pi/6*rho*dVec.^3;
    kern=coag_kernel(mp,mp,dVec,dVec,T,visc);
    
    summa = w'*kern'*w;

    if abs(c)<1e-3
        out=-pi*summa;
    else
        out = -pi*summa*c^2/(exp(c)-1)^2;
    end


end
