function out = coagS0lossTo0(bins,alpha,d1,d2,rho,T,visc)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    
   
   
    if d2/d1-1<1e-3
        out=-6.2832*d1^2*coag_kernel(pi/6*rho*d1^3,pi/6*rho*d2^3,d1,d2,T,visc);
        return;
    end
    
    
    dVec=logspace(lg(d1),lg(d2),bins);
    
    mp=pi/6*rho*dVec.^3;
    kern=coag_kernel(mp,mp,dVec,dVec,T,visc);
    
    [~,dp1]=meshgrid(dVec);
    kerroin=-dp1.^2;
    kern=kern.*kerroin;
    
    nj=n_j_powerlaw(1,alpha,d1,d2,dVec);
    kern(1,:)=kern(1,:)/2;
    kern(end,:)=kern(end,:)/2;
    kern(:,1)=kern(:,1)/2;
    kern(:,end)=kern(:,end)/2;
    
    out = nj*kern'*nj'*(ln(dVec(2)/dVec(1)))^2*2*pi;

%     for i=1:4
%         dp1=dVec(i);
%         
%         mp1=pi/6*rho*dp1^3;
% 
%         
%         summa2=0;
%         for j=1:4
%             dp2=dVec(j);
%             
%             mp2=pi/6*rho*dp2^3;
%   
%            
%             kern=coag_kernel(mp1,mp2,dp1,dp2,T,visc);
%             
%             summa2=summa2-dp1^2*kern*GaussOlinWeights4(j,c);
%             
%         end
%         
%         summa=summa+summa2*GaussOlinWeights4(i,c);
%     end




end
