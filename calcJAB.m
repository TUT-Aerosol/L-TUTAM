function fVec = calcJAB( aVec )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fVec=zeros(2,2);

a=aVec(1);
d=aVec(2);

if abs(a)<1e-6 || abs(a+2)<1e-6 || abs(a+3)<1e-6  
    a=a+1e-6;
end

if abs(d-1)<1e-6
    d=d+1e-6;
end

a2=a+2;
a3=a+3;
da1=d^(a+1);
da2=d^a2;
da3=d^a3;
lnd=log(d);
da=d^a;



% dA/da

A=((da3*(1+a*lnd)-1)*(da*a3-a3)-(da*(1+a3*lnd)-1)*(a*(da3-1)))/(da*a3-a3)^2;

    
fVec(1,1)=A/3/(a/a3*(1-da3)/(1-da))^(2/3);

% dA/dd

A=a/a3*(a3*da2*(da-1)-a*d^(a-1)*(da3-1))/(da-1)^2;
    
fVec(1,2)=A/3/(a/a3*(1-da3)/(1-da))^(2/3);

% dB/da
B=((da3*(1+a2*lnd)-1)*(da2*a3-a3)-(da2*(1+a3*lnd)-1)*(a2*(da3-1)))/(da2*a3-a3)^2;

    
fVec(2,1)=B;

% dB/dd
B=a2/a3*(a3*da2*(da2-1)-a2*da1*(da3-1))/(da2-1)^2;
    
fVec(2,2)=B;

if ~isreal(fVec)
    disp(aVec)
    error('imag JAB')
end

end

