function out=condS0(GR,a,d)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    if abs(d-1)<1e-3
        out=2*GR;
    elseif abs(a+1)<1e-3
        out=2*GR*ln(d)/(1-1/d);
    elseif abs(a+2)<1e-3
        out=2*GR*(d-1)/ln(d);
    else
        out=2*GR*(a+2)/(a+1)*(1-d^(-a-1))/(1-d^(-a-2));
    end

end

