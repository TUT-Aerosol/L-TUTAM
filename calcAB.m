function fVec = calcAB( aVec )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fVec=zeros(2,1);

a=aVec(1);
d=1/aVec(2);

if abs(a)<1e-6 || abs(a+2)<1e-6 || abs(a+3)<1e-6  
    a=a+1e-6;
end

if abs(d-1)<1e-6
    d=d+1e-6;
end



A=a/(a+3)*(1-d^(-3-a))/(1-d^(-a));

    
fVec(1)=A^(1/3);



B=(a+2)/(a+3)*(1-d^(-3-a))/(1-d^(-2-a));

    
fVec(2)=B;

if ~isreal(fVec)
    error('imag AB')
end

end

