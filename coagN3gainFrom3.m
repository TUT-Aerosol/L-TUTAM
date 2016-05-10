function out = coagN3gainFrom3(bins,rho,T,visc,cmd,ln2s)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    
    summa = 0;
    
    
    sigma = exp(sqrt(ln2s));
    dlogVec=logspace(lg(cmd*sigma^(-3)),lg(cmd*sigma^(3)),bins);
    dlndp=ln(dlogVec(2)/dlogVec(1));

    for i=1:bins
        dp1=dlogVec(i);
        
        mp1=pi/6*rho*dp1^3;
        nj1=n_j(dp1,1,sigma,cmd);
        
        if i==1 || i==bins
            nj1=nj1/2;
        end   
        
        summa2=0;
        for j=1:bins
            dp2=dlogVec(j);
            
            mp2=pi/6*rho*dp2^3;
            nj2=n_j(dp2,1,sigma,cmd);
            
            if j==1 || j==bins
                nj2=nj2/2;
            end    
           
            kern=coag_kernel(mp1,mp2,dp1,dp2,T,visc);
            
            summa2=summa2+kern*nj2;
            
        end
        
        summa=summa+summa2*nj1;
    end

    out = pi*summa*dlndp^2;

end
