function  nj  = n_j( dp,ntot,sigma,cmd )
%N_J dN/DlnDp

nj=zeros(size(dp));
for i=1:length(ntot)
    ala = sqrt(2*pi)*log(sigma(i));
	expoyla = (log(dp/cmd(i))).^2;
	expoala = 2*(log(sigma(i)))^2;
	expo = -1.0*expoyla/expoala;
	nj = nj+ntot(i)/ala*exp(expo);
end
    
        
end

