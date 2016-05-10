function [gmdVec,gsdVec] = laskeGmdGsdNumeerisesti(N0Vec,alpha0Vec,d1,dd0Vec,N3Vec,cmdVec,ln2sVec)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

gmdVec = zeros(size(N0Vec));
gsdVec = zeros(size(N0Vec));

for i=1:numel(N0Vec)
    NPL=N0Vec(i);
    NLN=N3Vec(i);
    N=NPL+NLN;
    dd=dd0Vec(i);
    d2=d1/dd;
    alpha=alpha0Vec(i);
    cmd=cmdVec(i);
    ln2s=ln2sVec(i);
    
    if NPL<1e-4
        gmdVec=cmdVec;
        
        dg=gmdVec(i);
        dp=logspace(lg(d1),-6,1000);
        nj=n_j(dp,NLN,exp(sqrt(ln2s)),cmd)+n_j_powerlaw(NPL,alpha,d1,d2,dp);
        muuttuja=(ln(dp./dg)).^2.*nj*ln(dp(2)/dp(1));
        ln2gsd=sum(muuttuja)/N;
        gsdVec(i)=exp(sqrt(ln2gsd));
    else
        
    
        if N>1e-4
            lndg=NPL/N*(ln(d2)-dd^alpha*ln(d1))/(1-dd^alpha)-NPL/N/alpha+NLN/(NPL+NLN)*ln(cmd);
            dg=exp(lndg);

            dp=logspace(lg(d1),-6,1000);
            nj=n_j(dp,NLN,exp(sqrt(ln2s)),cmd)+n_j_powerlaw(NPL,alpha,d1,d2,dp);
            muuttuja=(ln(dp./dg)).^2.*nj*ln(dp(2)/dp(1));
            ln2gsd=sum(muuttuja)/N;

        else
            dg=1e-9;
            ln2gsd=0;
        end

            gmdVec(i)=dg;
            gsdVec(i)=exp(sqrt(ln2gsd));
    end
    
end


end

