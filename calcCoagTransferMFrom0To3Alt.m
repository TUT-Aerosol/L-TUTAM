function out = calcCoagTransferMFrom0To3Alt(bins,alpha,d1,d2,rho,T,visc)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    
    if  d2/d1 > 3
        out = calcCoagTransferMFrom0To3(bins,alpha,d1,d2,rho,T,visc);
        return;
    end

    summa = 0;
    
    
    if d2/d1-1<1e-3
        out=d2^3*2*pi*coag_kernel(pi/6*rho*d1^3,pi/6*rho*d2^3,d1,d2,T,visc);
        return;
    end
    
    c=alpha*ln(d2/d1);
    
    daVec=logspace(lg(d1),lg(d2),4);

    



    for i=1:4
        dp1=daVec(i);
        
        mp1=pi/6*rho*dp1^3;
        
        dp_alku_b = alaraja(d2^3-dp1^3,d1^3)^(1/3);
        dbVec=logspace(lg(dp_alku_b),lg(d2),4);

        c2=alpha*ln(d2/dp_alku_b);
        
        w=GaussOlinWeights4((1:4)',c2);
        mp2=pi/6*rho*dbVec.^3;
        kern=coag_kernel(mp1,mp2,dp1,dbVec,T,visc);
        kerroin = (dp1^3+dbVec.^3);
        kern=kern.*kerroin;
        summa2 = kern*w;
        
        
%         summa2=0;
%         for j=1:4
%             dp2=dbVec(j);
%             
%             mp2=pi/6*rho*dp2^3;
% 
%             kern=coag_kernel(mp1,mp2,dp1,dp2,T,visc);
%             
%             summa2=summa2+((dp1^3+dp2^3)^(2/3))*kern*GaussOlinWeights4(j,c2);
%             
%         end
        
        if abs(c2)<1e-3
            summa=summa+summa2*GaussOlinWeights4(i,c)*ln(d2/dp_alku_b)/ln(d2/d1);
        else
            summa=summa+summa2*GaussOlinWeights4(i,c)*c2/(exp(c2)*(1-(d1/d2)^alpha));
        end
        
    end

    if abs(c)<1e-3
        out=pi*summa;
    else
        out = pi*summa*c/(exp(c)-1);
    end
    

    

end
