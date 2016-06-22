function [alphaLaskettu,dLaskettu] = laskeAlphaDTaulukoista(Aoikea,Boikea,p,alphaInit,dInit)

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
    
    if p.PLEquations == 1
        alphaLaskettu=p.alphaTaulukko(Aoikea,Boikea);
        dLaskettu=p.dTaulukko(Aoikea,Boikea);
        
    else
        [alphaLaskettu, dLaskettu] = LMAEquations(Aoikea,Boikea,alphaInit,dInit);
    end

    if ~isreal(alphaLaskettu)
        error('imag alpha')
    end
    
    if ~isreal(dLaskettu)
        error('imag d')
    end    
end