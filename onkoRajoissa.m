function onRajoissa = onkoRajoissa( aVec,alphaMin,alphaMax,dMin,dMax )
alpha=aVec(1);
d=aVec(2);

onRajoissa=1;

if alpha<alphaMin
    onRajoissa=0;
    return;
end

if alpha>alphaMax
    onRajoissa=0;
    return;
end

if d<dMin
    onRajoissa=0;
    return;
end

if d>dMax
    onRajoissa=0;
    return;
end


end

