function out = lossM0(a,d,n)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    
    if abs(d-1)<1e-3
        out = 1;
    else
        if abs(a+3)<1e-3
            out = (1-d^(-n))/(n*ln(d));
        elseif abs(a+n+3)<1e-3
            out = n*ln(d)/(d^n-1);
        else
            out = (a+3)/(a+n+3)*(1-d^(-n-3-a))/(1-d^(-a-3));
        end
    end


end

