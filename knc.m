function  knc  = knc( mp1,mp2,dp1,dp2,temp, visc,diff1,diff2 )
    
    velo1 = velo(mp1,temp);
    velo2 = velo(mp2,temp);
    
    knc = 6.0*(diff1+diff2)./((dp1+dp2).*sqrt(velo1^2+velo2.^2));
    
end

