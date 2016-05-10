function [alphaLaskettu,dLaskettu] = laskeAlphaDTaulukoista(Aoikea,Boikea,p)

    if Aoikea < 1
        alphaLaskettu=1;
        dLaskettu=1;
        return
    elseif Boikea < 1
        alphaLaskettu=1;
        dLaskettu=1;
        return
    elseif isinf(Aoikea)
        alphaLaskettu=1;
        dLaskettu=1;
        return
    elseif isinf(Boikea)
        alphaLaskettu=1;
        dLaskettu=1;
        return
    elseif isnan(Aoikea)
        alphaLaskettu=1;
        dLaskettu=1;
        return
    elseif isnan(Boikea)
        alphaLaskettu=1;
        dLaskettu=1;
        return
    end
    
    alphaLaskettu=p.alphaTaulukko(Aoikea,Boikea);
    dLaskettu=p.dTaulukko(Aoikea,Boikea);

end