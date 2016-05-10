function [ diff_p ] = diff_p( dp,temp, visc )

    k = 1.381e-23;
    cc = cunningham(dp,temp);
    
    diff_p = k.*temp.*cc./(3*pi*visc.*dp);
    
end

