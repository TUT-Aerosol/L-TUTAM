function h=plot_lognormal_powerlaw(N0,alpha0,d1,d2,ntot,cmd,sigma)
    
    dp=logspace(-9,-6,1000);

    nj_log=n_j(dp,ntot,sigma,cmd);
    nj_powerlaw=n_j_powerlaw(N0,alpha0,d1,d2,dp);
    
    

%     h=semilogx(dp*1e9,nj_log*ln(10),'-b',dp*1e9,nj_powerlaw*ln(10),'-g','linewidth',2); % erikseen
    h=semilogx(dp*1e9,(nj_powerlaw+nj_log)*ln(10),'-b','linewidth',2); % summana
    xlabel('dp (nm)')
    ylabel('dN/dlogdp (cm^{-3})')
    

end

