function  beta  = coag_kernel_dp( dp1,dp2,temp,rhoP )

    if nargin < 4
        rhoP=1000;
        
        if nargin < 3
            temp=300;
            
            if nargin < 2
                dp2=dp1;
            end
        end
    end
    
    if length(dp1)>1
        fprintf('coag_kernel_dp error: The first argument MUST be scalar!! \n')
        return
    end
    
    visc = viscosity(temp);
    
    dp1=max(dp1,1e-10);
    mp1=0.523598775598299*rhoP*dp1^3;
    diff1 = diff_p(dp1,temp,visc);

    dp2=max(dp2,1e-10);
    mp2=0.523598775598299*rhoP.*dp2.^3;
    diff2 = diff_p(dp2,temp,visc);
    trans = trans_dahneke(mp1,mp2,dp1,dp2,temp,visc,diff1,diff2);
    beta=trans.*(dp1+dp2).*(diff1+diff2);
        

end

