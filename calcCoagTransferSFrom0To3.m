function out = calcCoagTransferSFrom0To3(bins,alpha,d1,d2,rho,T,visc)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    if d2/d1-1<1e-3
        out=d2^2*4.986967483164005*coag_kernel(pi/6*rho*d1^3,pi/6*rho*d2^3,d1,d2,T,visc);
        return;
    end

    summa = 0;
    
    
    daVec=logspace(lg(d1),lg(d2),bins);
    nja=n_j_powerlaw(1,alpha,d1,d2,daVec);
    nja([1 end])=nja([1 end])/2;

    for i=1:bins
        dp1=daVec(i);
        
        mp1=pi/6*rho*dp1^3;
        
        dp_alku_b = alaraja(d2^3-dp1^3,d1^3)^(1/3);
        dbVec=logspace(lg(dp_alku_b),lg(d2),bins);
        njb=n_j_powerlaw(1,alpha,d1,d2,dbVec);
        njb([1 end])=njb([1 end])/2;

        mp2=pi/6*rho*dbVec.^3;
        kern=coag_kernel(mp1,mp2,dp1,dbVec,T,visc);
        kerroin = (dp1^3+dbVec.^3).^(2/3);
        kern=kern.*kerroin;
        summa2 = kern*njb';
        
%         summa2=0;
%         for j=1:4
%             dp2=dbVec(j);
%             
%             mp2=pi/6*rho*dp2^3;
% 
%             kern=coag_kernel(mp1,mp2,dp1,dp2,T,visc);
%             
%             summa2=summa2+kern*GaussOlinWeights4(j,c2);
%             
%         end
        

        summa=summa+summa2*nja(i)*ln(dbVec(2)/dbVec(1));


        
    end

    out=pi*summa*ln(daVec(2)/daVec(1));
    

    

end
