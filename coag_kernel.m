function  beta  = coag_kernel( mp1Vec,mp2,dp1Vec,dp2,temp, visc )

%     beta=2.2e-16;return;

    beta=zeros(length(dp1Vec),length(dp2));
    dp2=max(dp2,1e-10);
    diff2 = diff_p(dp2,temp,visc);
    
    
    for i=1:length(dp1Vec)
        
        dp1=max(dp1Vec(i),1e-10);
        mp1=mp1Vec(i);

        diff1 = diff_p(dp1,temp,visc);


        trans = trans_dahneke(mp1,mp2,dp1,dp2,temp,visc,diff1,diff2);

        beta(i,:)=trans.*(dp1+dp2).*(diff1+diff2);
    end
    
    
end

