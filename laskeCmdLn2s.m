function [cmd,ln2s]=laskeCmdLn2s(N,S,M)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
    cmd=zeros(size(N));
    ln2s=cmd;
    
    for i=1:numel(N)

        if N(i)<1e-8
            cmd(i)=1.1e-9;
            ln2s(i)=log(1.01)^2;

        elseif S(i)<1e-30
            cmd(i)=1.1e-9;
            ln2s(i)=log(1.01)^2; 

        elseif M(i)<1e-35
            cmd(i)=1.1e-9;
            ln2s(i)=log(1.01)^2;
                        
        else

            ln2s(i) = 2/3*log(M(i)) + 1/3*log(N(i)) - log(S(i));
            ln2s(i) = rajat(log(1.01)^2,ln2s(i),log(2)^2);

            cmd(i) = alaraja((M(i)/N(i))^(1/3)*exp(-1.5*ln2s(i)),1.1e-9);
        end
        
        if isnan(cmd(i))
            error('nan cmd');
        end
        
        if isinf(cmd(i))
            error('inf cmd');
        end
    
    end

end

