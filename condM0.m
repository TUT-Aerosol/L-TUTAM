function out=condM0(GR,a,d)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    if abs(d-1)<1e-3
        out=3*GR;
    elseif abs(a+2)<1e-3
        out=3*GR*ln(d)/(1-1/d);
    elseif abs(a+3)<1e-3
        out=3*GR*(d-1)/ln(d);
    else
        out=3*GR*(a+3)/(a+2)*(1-d^(-a-2))/(1-d^(-a-3));
    end

end

