function out = lossN0(a,d,n)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    
    if abs(d-1)<1e-3
        out = 1;
    else
        if abs(a)<1e-3
            out = (1-d^(-n))/(n*ln(d));
        elseif abs(a+n)<1e-3
            out = n*ln(d)/(d^n-1);
        else
            out = a/(a+n)*(1-d^(-n-a))/(1-d^(-a));
        end
    end


end

