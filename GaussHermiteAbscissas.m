function out=GaussHermiteAbscissas(n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
out = zeros(size(n));

for i=1:numel(n)
    if n(i)==1
        out(i)=-2.020182870456086;
    elseif n(i)==2
        out(i)=-0.958572464613819;
    elseif n(i)==3
        out(i)=0;
    elseif n(i)==4
        out(i)=0.958572464613819;
    else
        out(i)=2.020182870456086;
    end
end

end

