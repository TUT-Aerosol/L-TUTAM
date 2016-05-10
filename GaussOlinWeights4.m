function out=GaussOlinWeights4(i,c)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    c2=c^2;
    c3=c^3;
    c4=c^4;
    ec=exp(c);

    out = zeros(size(i));
    for ind=1:numel(i)
        if abs(c)<1e-3
            if i(ind)==1 || i(ind)==4
                out(ind)=1/8;
            else
                out(ind)=3/8;
            end
            
        else
        

            if i(ind)==1
                out(ind)=(27*ec-27)/c4+(-9*ec-18)/c3+(ec-5.5)/c2-1/c;
            elseif i(ind)==2
                out(ind)=(-81*ec+81)/c4+(36*ec+45)/c3+(-4.5*ec+9)/c2;
            elseif i(ind)==3
                out(ind)=(81*ec-81)/c4-(45*ec+36)/c3+(9*ec-4.5)/c2;
            else
                out(ind)=(-27*ec+27)/c4+(18*ec+9)/c3+(-5.5*ec+1)/c2+ec/c;
            end
        end
    end

end

