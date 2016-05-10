function out = coagM0lossTo0(bins,alpha,d1,d2,rho,T,visc)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    if d2/d1-1<1e-3
        out=-6.2832*(d1^3)*coag_kernel(pi/6*rho*d1^3,pi/6*rho*d2^3,d1,d2,T,visc);
        return;
    end
    
    dVec=logspace(lg(d1),lg(d2),bins);
    
    mp=pi/6*rho*dVec.^3;
    kern=coag_kernel(mp,mp,dVec,dVec,T,visc);
    
    [~,dp1]=meshgrid(dVec);
    kerroin=-dp1.^3;
    kern=kern.*kerroin;
    
       nj=n_j_powerlaw(1,alpha,d1,d2,dVec);
    kern(1,:)=kern(1,:)/2;
    kern(end,:)=kern(end,:)/2;
    kern(:,1)=kern(:,1)/2;
    kern(:,end)=kern(:,end)/2;
    
    out = nj*kern'*nj'*(ln(dVec(2)/dVec(1)))^2*2*pi;

end
