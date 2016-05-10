function out = coagM3lossTo0(bins,alpha,d1,d2,rho,T,visc,cmd,ln2s)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    
    if d2/d1-1<1e-3
        out=-2*pi*cmd^3*coag_kernel(pi/6*rho*d1^3,pi/6*rho*cmd^3,d1,cmd,T,visc);
        return;
    end
    
    dVec=logspace(lg(d1),lg(d2),bins);
    mp1=pi/6*rho*dVec.^3;
    
    x2=GaussHermiteAbscissas((1:5)');
    w2=GaussHermiteWeights((1:5)');
%     [x2,w2]=GaussHermite(bins);
    d2Vec=cmd*exp(x2*sqrt(2*ln2s));
    mp2=pi/6*rho*d2Vec.^3;
    
    kern=coag_kernel(mp1,mp2,dVec,d2Vec,T,visc);
    
    [dp2,~]=meshgrid(d2Vec,dVec);
    kerroin=dp2.^3;
    kern=kern.*kerroin;
    
    nj=n_j_powerlaw(1,alpha,d1,d2,dVec);
    kern(1,:)=kern(1,:)/2;
    kern(end,:)=kern(end,:)/2;
%     kern(:,1)=kern(:,1)/2;
%     kern(:,end)=kern(:,end)/2;
    
%     out = -nj*kern*w2*(ln(dVec(2)/dVec(1)))*2*sqrt(pi);
    out = -w2'*kern'*nj'*(ln(dVec(2)/dVec(1)))*2*sqrt(pi);


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
% 
%     end

    
end
