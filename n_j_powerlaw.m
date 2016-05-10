function  nj  = n_j_powerlaw( N,alpha,d1,d2,dp )
%N_J dN/DlnDp

nj=zeros(size(dp));

if (d2/d1-1)<1e-3
    return
else

    for i=1:numel(dp)
        if dp(i)<d1-1e-12
            nj(i)=0;
        elseif dp(i)>d2+1e-12
            nj(i)=0;
        else
            if abs(alpha)<1e-3
                nj(i)=N/ln(d2/d1);
            else
                nj(i)=N*alpha./((d2./dp(i)).^alpha-(d1./dp(i)).^alpha);
            end
        end
    end
end
        
    
end

