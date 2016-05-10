function out=GaussHermiteWeights(n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
out=zeros(size(n));

for i=1:numel(n)
    if n(i)==1 || n(i)==5
        out(i)=0.019953242059;
    elseif n(i)==2 || n(i)==4
        out(i)=0.393619323152;
    else
        out(i)=0.945308720483;
    end
end

end

