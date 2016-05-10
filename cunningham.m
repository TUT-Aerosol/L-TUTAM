function [ cc ] = cunningham( dp,temp )

    lambda = vapaa_matka(temp);
    
    expo = -0.39.*dp./lambda;
    kerroin = 2.34+1.05.*exp(expo);
    
    cc = 1+lambda./dp.*kerroin;
    
end

