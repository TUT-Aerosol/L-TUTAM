function out = coagN3lossTo0Alt(bins,alpha,d1,d2,rho,T,visc,cmd,ln2s)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    
    if  d2/d1 > 3
        out = coagN3lossTo0(bins,alpha,d1,d2,rho,T,visc,cmd,ln2s);
        return;
    end
   
    if d2/d1-1<1e-3
        out=-2*pi*coag_kernel(pi/6*rho*d1^3,pi/6*rho*cmd^3,d1,cmd,T,visc);
        return;
    end
    
    c=alpha*ln(d2/d1);
    w=GaussOlinWeights4((1:4)',c);
    
    dVec=logspace(lg(d1),lg(d2),4);
    mp1=pi/6*rho*dVec.^3;
    
    x2=GaussHermiteAbscissas((1:5)');
    w2=GaussHermiteWeights((1:5)');
    d2Vec=cmd*exp(x2*sqrt(2*ln2s));
    mp2=pi/6*rho*d2Vec.^3;
    
    kern=coag_kernel(mp1,mp2,dVec,d2Vec,T,visc);
    
    summa = w2'*kern'*w;
        

%     for i=1:4
%         dp1=dVec(i);
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
%             summa2=summa2+kern*GaussHermiteWeights(j);
%             
%         end
%         
%         summa=summa+summa2*GaussOlinWeights4(i,c);
%     end

    if abs(c)<1e-3
        out=-2*sqrt(pi)*summa;
    else
        out =-2*sqrt(pi)*summa*c/(exp(c)-1);
    end
    
    
end
