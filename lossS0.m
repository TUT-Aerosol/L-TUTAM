function out = lossS0(a,d,n)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    
    if abs(d-1)<1e-3
        out = 1;
    else
        if abs(a+2)<1e-3
            out = (1-d^(-n))/(n*ln(d));
        elseif abs(a+n+2)<1e-3
            out = n*ln(d)/(d^n-1);
        else
            out = (a+2)/(a+2+n)*(1-d^(-n-2-a))/(1-d^(-a-2));
        end
    end


end

