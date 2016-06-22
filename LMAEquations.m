function [ alpha, d ] = LMAEquations( A,B,alphaInit,dInit )
    A=A^(1/3);

    if A>B || A<1 || B<1
        d=1;
        alpha=1;
        return;
    end
    
    if B>87 || A>87
        d=1/100;
        alpha=5;
        return;
    end
    
    yOikea=[A;B];

    aInit=[alphaInit;1/dInit];
    alphaMin=-10000;
    alphaMax=10000;
    dMin=0;
    dMax=1000;

    lambda=0.001;
    lambdaDown=10;
    lambdaUp=5;
    convCriterion=0.001;

    aNyt=aInit;

    yInit=calcAB(aInit);
    yNyt=yInit;
    i=1;
    iUp=0;

    iterMax=40;


    while max(abs(yNyt./yOikea-1))>convCriterion



        J=calcJAB(aNyt);
        JT=J';
        JTJ=JT*J;
        r=yOikea-calcAB(aNyt);
        costGradient=JT*r;
        g=JTJ+lambda*eye(2);
        costNyt=sum(r.^2);

        delta=g\costGradient;

        aUusi=aNyt+delta;

        rUusi=yOikea-calcAB(aUusi);
        costUusi=sum(rUusi.^2);

        if costUusi<costNyt && onkoRajoissa(aUusi,alphaMin,alphaMax,dMin,dMax) && ~any(any(isnan(calcJAB(aUusi))))
            aNyt=aUusi;
            lambda=lambda/lambdaDown;

            yNyt= calcAB(aNyt);
            

            i=i+1;


        else
            iUp=iUp+1;
            lambda=lambda*lambdaUp;
            if lambda>1e10
                break;
            end
        end

    end

    if max(abs(yNyt./yOikea-1))>convCriterion
        alpha=0;
        d=1;
        error eiConverg
    else

        alpha=max(min(5,aNyt(1)),-5);
        d=max(min(1,1/aNyt(2)),0.01);
    end


end

