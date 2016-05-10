function out = coagM3gainFrom3Alt(bins,rho,T,visc,cmd,ln2s)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    
%     summa = 0;
    
    w=GaussHermiteWeights((1:5)');
    x=GaussHermiteAbscissas((1:5)');
    dVec=cmd*exp(x*sqrt(2*ln2s));
    
    mp=pi/6*rho*dVec.^3;
    kern=coag_kernel(mp,mp,dVec,dVec,T,visc);
    
    [dp2,dp1]=meshgrid(dVec);
    kerroin=dp1.^3+dp2.^3;
    kern=kern.*kerroin;
    
    summa = w'*kern'*w;

%     for i=1:5
%         x1=GaussHermiteAbscissas(i);
%         dp1=cmd*exp(x1*sqrt(2*ln2s));
%         
%         mp1=pi/6*rho*dp1^3;
% 
%         
%         summa2=0;
%         for j=1:5
%             x2=GaussHermiteAbscissas(j);
%             dp2=cmd*exp(x2*sqrt(2*ln2s));
%             
%             mp2=pi/6*rho*dp2^3;
% 
%             kern=coag_kernel(mp1,mp2,dp1,dp2,T,visc);
%             
%             summa2=summa2+(dp1^3+dp2^3)*kern*GaussHermiteWeights(j);
%             
%         end
%         
%         summa=summa+summa2*GaussHermiteWeights(i);
%     end

    out = summa;

end
