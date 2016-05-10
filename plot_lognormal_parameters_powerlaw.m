function h=plot_lognormal_parameters_powerlaw(N0,alpha0,d1,d2,N3,cmd,ln2s) 


    h=plot_lognormal_powerlaw(N0,alpha0,d1,d2,N3,cmd,exp(sqrt(ln2s)));

end

