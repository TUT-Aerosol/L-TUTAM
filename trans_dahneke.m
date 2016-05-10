function  trans  = trans_dahneke_oikeasti( mp1,mp2,dp1,dp2,temp, visc,diff1,diff2 )

    kncee = 2/3*knc(mp1,mp2,dp1,dp2,temp,visc,diff1,diff2);
    
    trans = (1+kncee)./(1+2*kncee+2*kncee.^2);
    
end

