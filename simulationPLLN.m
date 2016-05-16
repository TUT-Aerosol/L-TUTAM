function out =  simulationPLLN(p)
%SIMULATIONPLLN simulates the PLLN model using input parameter set 'p'    

    % 2nd and 3rd moments of a newly formed particle
    p.sCluster = p.dCluster^2;
    p.mCluster = p.dCluster^3;

    % Distributions are sometimes indexed as numbers that are:
    % 0: PL
    % 3: LN
    switch(p.model)
        case 'PLLN'
            p.modes=30;
        case 'LN'
            p.modes=3;
    end

    % Total simulation time
    p.totalTime = p.timeVec(end)-p.timeVec(1);
    
    % If 'p.initialMomentVec' is scalar 0, make a vector
    if length(p.initialMomentVec)==1
        if p.initialMomentVec==0
            p.initialMomentVec = zeros(1,6);
        end
    end
    
    % Calculating coagulation sink factor
    p.coagSinkFactor = coag_kernel_dp(5e-9,p.coagSinkCMD,p.T,p.rho)*2*pi*5e-9^(-p.coagSinkExponent);
    
    % Options to ODE solver
    options = odeset('RelTol',p.relativeTolerance,'nonnegative',1:6,'stats','off');
    
    if p.plotWaitbarDuringSim
        hWait = waitbar(0,'Simulating...');
    end
    
    % Solving ODE
    [t,Y] = feval(p.solverName,@modelFunc,p.timeVec,p.initialMomentVec,options,p);

    if p.plotWaitbarDuringSim
        close(hWait);
    end
    
    disp('Simulation performed')
    disp('Making output struct...')
    
    % Making the output struct 'out'
    % Parameters of the distribution need to be computed from the moment data
    out.t = t; 
    out.Y = Y;
    out.p = p;
    out.N_LN = Y(:,4);
    [cmdVec,ln2sVec] = laskeCmdLn2s(Y(:,4),Y(:,5),Y(:,6));
    
    out.CMD = cmdVec;
    out.sigma = exp(sqrt(ln2sVec));

    alpha0Vec=zeros(size(Y(:,1)));
    dd0Vec=alpha0Vec;
    for i=1:length(Y(:,1))
        [alpha0Vec(i),dd0Vec(i)]=laskeAlphaDTaulukoista(Y(i,3)/Y(i,1)/p.mCluster,Y(i,3)/Y(i,2)/p.dCluster,p);
    end
    
    out.N_PL = Y(:,1);
    out.alpha = alpha0Vec;
    out.D2 = p.dCluster./dd0Vec;
    out.N = Y(:,1)+Y(:,4);
    out.M_2 = Y(:,2)+Y(:,5);
    out.M_3 = Y(:,3)+Y(:,6);
    
    % Computing GMD and GSD of the total distr.
    [gmdVec,gsdVec] = laskeGmdGsdNumeerisesti(Y(:,1),alpha0Vec,p.dCluster,dd0Vec,Y(:,4),cmdVec,ln2sVec);
    out.GMD = gmdVec;
    out.GSD = gsdVec;

    % Removing a large interpolation table from the output struct
    out.p=rmfield(out.p,'alphaTaulukko');
    out.p=rmfield(out.p,'dTaulukko');

end


function dy = modelFunc(t,y,param)
%MODELFUNC is the internal function containing the PLLN model equations

    if any(isnan(y))
        disp(y)
        error('NaN in y')
    end
    
    if param.plotWaitbarDuringSim
        waitbar(t/param.totalTime);
    end
    
    % Compute the parameters of the distributions
    [cmd,ln2s]=laskeCmdLn2s(y(4),y(5),y(6));

    d1_0=param.dCluster;
    m1_0=d1_0^3;
        
    A=y(3)/y(1)/m1_0;
    B=y(3)/y(2)/d1_0;
    [alpha0,dd0] = laskeAlphaDTaulukoista(A,B,param);
    
    d2_0=param.dCluster/dd0;
    
    % Find J and GR from the matrices
    param.J=param.JMatrix(2,find(t<=param.JMatrix(1,:),1));
    param.GR=param.GRMatrix(2,find(t<=param.GRMatrix(1,:),1));
    
    % Plotting the distr.
    if param.plotDistrDuringSim
        figure(1)
        clf
        plot_lognormal_parameters_powerlaw(y(1),alpha0,d1_0,d2_0,y(4),cmd,ln2s);

        ylim([max(param.JMatrix(2,:))/100 max(param.JMatrix(2,:))*10]*param.totalTime)
        xlim([1 100])
        set(gca,'yscale','log')
        title(strcat('t= ',num2str(t,'%1.3f'),' s'))
        xlabel('Dp (nm)')
        ylabel('dN/dlogDp (cm^{-3})')

        drawnow
    end

    
    dy = zeros(size(y));
    % 1 mode0 N
    % 2 mode0 S
    % 3 mode0 M

    % 4 mode3 N
    % 5 mode3 S
    % 6 mode3 M

    if param.modes == 0 || param.modes == 30
        % new particle formation to PL
        dy(1) = dy(1) + param.J;
        dy(2) = dy(2) + param.J*param.sCluster;
        dy(3) = dy(3) + param.J*param.mCluster;
        
        % condensation in PL
        if param.sizeDependentGR == 0
            dy(2) = dy(2) + condS0(param.GR,alpha0,dd0)*y(2)/d1_0;
            dy(3) = dy(3) + condM0(param.GR,alpha0,dd0)*y(3)/d1_0;
        else
            if param.numericCondensation == 1 || alpha0 < 0.5
                dy(2) = dy(2) + condS0SizeDependentNumeric(param.condBins,alpha0,d1_0,d2_0,param.dsa,param.diffsa,param.CsaInf,param.CsaSat,param.rh,param.T,param.visc)*y(1);
                dy(3) = dy(3) + condM0SizeDependentNumeric(param.condBins,alpha0,d1_0,d2_0,param.dsa,param.diffsa,param.CsaInf,param.CsaSat,param.rh,param.T,param.visc)*y(1);
            else
                dy(2) = dy(2) + condS0SizeDependent(alpha0,d1_0,d2_0,param.dsa,param.diffsa,param.CsaInf,param.CsaSat,param.rh,param.T,param.visc)*y(1);
                dy(3) = dy(3) + condM0SizeDependent(alpha0,d1_0,d2_0,param.dsa,param.diffsa,param.CsaInf,param.CsaSat,param.rh,param.T,param.visc)*y(1);
            end
        end 
    end
    
    if param.modes == 30
        % condensational transfer PL -> LN
        if param.condensationalTransfer
            if param.sizeDependentGR == 0
                GRD2 = param.GR;
            else
                GRD2 = massGrowthRate(d2_0,param.dsa,param.T,param.visc,param.diffsa,param.CsaInf,param.CsaSat,param.rh)/(pi/2*rhoo(d2_0,param.T,param.rh)*d2_0^2);
            end
            n_j_powerlaw_d2_0 = n_j_powerlaw(y(1),alpha0,d1_0,d2_0,d2_0)*param.condensationalTransferFactor;
            dy(1) = dy(1) - GRD2/d2_0*n_j_powerlaw_d2_0;
            dy(2) = dy(2) - GRD2*d2_0*n_j_powerlaw_d2_0;
            dy(3) = dy(3) - GRD2*d2_0^2*n_j_powerlaw_d2_0;
            dy(4) = dy(4) + GRD2/d2_0*n_j_powerlaw_d2_0;
            dy(5) = dy(5) + GRD2*d2_0*n_j_powerlaw_d2_0;
            dy(6) = dy(6) + GRD2*d2_0^2*n_j_powerlaw_d2_0;
            
        end
    end
    
    if param.modes == 3
        % new particle formation to LN
        dy(4) = dy(4) + param.J;
        dy(5) = dy(5) + param.J*param.sCluster;
        dy(6) = dy(6) + param.J*param.mCluster;
    end
    
    if param.modes == 3 || param.modes == 30
        % condensation in LN
        if param.sizeDependentGR == 0
            dy(5) = dy(5) + 2*cmd^(-1)*exp(-3/2*ln2s)*param.GR*y(5);
            dy(6) = dy(6) + 3*cmd^(-1)*exp(-5/2*ln2s)*param.GR*y(6);
        else
            dy(5) = dy(5) + condS3SizeDependent(cmd,ln2s,param.dsa,param.diffsa,param.CsaInf,param.CsaSat,param.rh,param.T,param.visc)*y(4);
            dy(6) = dy(6) + condM3SizeDependent(cmd,ln2s,param.dsa,param.diffsa,param.CsaInf,param.CsaSat,param.rh,param.T,param.visc)*y(4);
        end
    end
  
    % losses
    if param.losses == 1
        dy(1) = dy(1) - param.lossesCoeff*d1_0^param.lossesExponent*lossN0(alpha0,dd0,param.lossesExponent)*y(1);
        dy(2) = dy(2) - param.lossesCoeff*d1_0^param.lossesExponent*lossS0(alpha0,dd0,param.lossesExponent)*y(2);
        dy(3) = dy(3) - param.lossesCoeff*d1_0^param.lossesExponent*lossM0(alpha0,dd0,param.lossesExponent)*y(3);
        dy(4) = dy(4) - param.lossesCoeff*cmd^param.lossesExponent*exp(param.lossesExponent^2/2*ln2s)*y(4);
        dy(5) = dy(5) - param.lossesCoeff*cmd^param.lossesExponent*exp((param.lossesExponent^2/2+param.lossesExponent*2)*ln2s)*y(5);
        dy(6) = dy(6) - param.lossesCoeff*cmd^param.lossesExponent*exp((param.lossesExponent^2/2+param.lossesExponent*3)*ln2s)*y(6);
    end

    % coagulation and coagulational transfer
    if param.coag == 1
        if param.numericCoagulation == 0
            
            % gaussian integration for PL, but numeric when d > 3,
            % gaussian for LN
            if param.modes == 0
                dy(1) = dy(1) + coagN0Alt(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc)*y(1)^2*1e6;
                dy(2) = dy(2) + coagS0Alt(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc)*y(1)^2*1e6;

            elseif param.modes == 30

                if param.coagulationalTransfer == 1
                    coagTransferNFrom0To3 = calcCoagTransferNFrom0To3Alt(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc)*y(1)^2*1e6;
                    coagTransferSFrom0To3 = calcCoagTransferSFrom0To3Alt(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc)*y(1)^2*1e6;
                    coagTransferMFrom0To3 = calcCoagTransferMFrom0To3Alt(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc)*y(1)^2*1e6;
                else
                    coagTransferNFrom0To3 = 0;
                    coagTransferSFrom0To3 = 0;
                    coagTransferMFrom0To3 = 0;
                end



                dy(1) = dy(1) + (coagN0gainFrom0Alt(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc)+ ...
                    coagN0lossTo0Alt(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc))*y(1)^2*1e6+ ...
                    coagN0lossTo3Alt(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc,cmd,ln2s)*y(1)*y(4)*1e6- ...
                    coagTransferNFrom0To3;

                dy(2) = dy(2) + (coagS0gainFrom0Alt(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc)+ ...
                    coagS0lossTo0Alt(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc))*y(1)^2*1e6+ ...
                    coagS0lossTo3Alt(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc,cmd,ln2s)*y(1)*y(4)*1e6- ...
                    coagTransferSFrom0To3;

                dy(3) = dy(3) + (coagM0gainFrom0Alt(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc)+ ...
                    coagM0lossTo0Alt(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc))*y(1)^2*1e6+ ...
                    coagM0lossTo3Alt(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc,cmd,ln2s)*y(1)*y(4)*1e6- ...
                    coagTransferMFrom0To3;

                dy(4) = dy(4) + coagN3gainFrom0Alt(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc,cmd,ln2s)*y(1)*y(4)*1e6 + ...
                     coagN3gainFrom3Alt(param.binsInCoagulation,param.rho,param.T,param.visc,cmd,ln2s)*y(4)^2*1e6 + ...
                     coagN3lossTo0Alt(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc,cmd,ln2s)*y(1)*y(4)*1e6 + ...
                     coagN3lossTo3Alt(param.binsInCoagulation,param.rho,param.T,param.visc,cmd,ln2s)*y(4)^2*1e6+ ...
                     coagTransferNFrom0To3;

                dy(5) = dy(5) + coagS3gainFrom0Alt(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc,cmd,ln2s)*y(1)*y(4)*1e6 + ...
                    coagS3gainFrom3Alt(param.binsInCoagulation,param.rho,param.T,param.visc,cmd,ln2s)*y(4)^2*1e6 + ...
                    coagS3lossTo0Alt(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc,cmd,ln2s)*y(1)*y(4)*1e6 + ...
                    coagS3lossTo3Alt(param.binsInCoagulation,param.rho,param.T,param.visc,cmd,ln2s)*y(4)^2*1e6+ ...
                    coagTransferSFrom0To3;

                dy(6) = dy(6) + coagM3gainFrom0Alt(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc,cmd,ln2s)*y(1)*y(4)*1e6 + ...
                    coagM3gainFrom3Alt(param.binsInCoagulation,param.rho,param.T,param.visc,cmd,ln2s)*y(4)^2*1e6 + ...
                    coagM3lossTo0Alt(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc,cmd,ln2s)*y(1)*y(4)*1e6 + ...
                    coagM3lossTo3Alt(param.binsInCoagulation,param.rho,param.T,param.visc,cmd,ln2s)*y(4)^2*1e6+ ...
                    coagTransferMFrom0To3;

            elseif param.modes == 3
                dy(4) = dy(4) + coagN3Alt(param.binsInCoagulation,cmd,ln2s,param.rho,param.T,param.visc)*y(4)^2*1e6;
                dy(5) = dy(5) + coagS3Alt(param.binsInCoagulation,cmd,ln2s,param.rho,param.T,param.visc)*y(4)^2*1e6;
            end
            
      % numeric integration for PL,
      % gaussian for LN     
        else
            if param.modes == 0
                dy(1) = dy(1) + coagN0(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc)*y(1)^2*1e6;
                dy(2) = dy(2) + coagS0(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc)*y(1)^2*1e6;

            elseif param.modes == 30

                if param.coagulationalTransfer == 1
                    coagTransferNFrom0To3 = calcCoagTransferNFrom0To3(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc)*y(1)^2*1e6;
                    coagTransferSFrom0To3 = calcCoagTransferSFrom0To3(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc)*y(1)^2*1e6;
                    coagTransferMFrom0To3 = calcCoagTransferMFrom0To3(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc)*y(1)^2*1e6;
                else
                    coagTransferNFrom0To3 = 0;
                    coagTransferSFrom0To3 = 0;
                    coagTransferMFrom0To3 = 0;
                end



                dy(1) = dy(1) + (coagN0gainFrom0(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc)+ ...
                    coagN0lossTo0(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc))*y(1)^2*1e6+ ...
                    coagN0lossTo3(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc,cmd,ln2s)*y(1)*y(4)*1e6- ...
                    coagTransferNFrom0To3;

                dy(2) = dy(2) + (coagS0gainFrom0(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc)+ ...
                    coagS0lossTo0(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc))*y(1)^2*1e6+ ...
                    coagS0lossTo3(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc,cmd,ln2s)*y(1)*y(4)*1e6- ...
                    coagTransferSFrom0To3;

                dy(3) = dy(3) + (coagM0gainFrom0(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc)+ ...
                    coagM0lossTo0(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc))*y(1)^2*1e6+ ...
                    coagM0lossTo3(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc,cmd,ln2s)*y(1)*y(4)*1e6- ...
                    coagTransferMFrom0To3;

                dy(4) = dy(4) + coagN3gainFrom0(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc,cmd,ln2s)*y(1)*y(4)*1e6 + ...
                     coagN3gainFrom3Alt(param.binsInCoagulation,param.rho,param.T,param.visc,cmd,ln2s)*y(4)^2*1e6 + ...
                     coagN3lossTo0(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc,cmd,ln2s)*y(1)*y(4)*1e6 + ...
                     coagN3lossTo3Alt(param.binsInCoagulation,param.rho,param.T,param.visc,cmd,ln2s)*y(4)^2*1e6+ ...
                     coagTransferNFrom0To3;

                dy(5) = dy(5) + coagS3gainFrom0(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc,cmd,ln2s)*y(1)*y(4)*1e6 + ...
                    coagS3gainFrom3Alt(param.binsInCoagulation,param.rho,param.T,param.visc,cmd,ln2s)*y(4)^2*1e6 + ...
                    coagS3lossTo0(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc,cmd,ln2s)*y(1)*y(4)*1e6 + ...
                    coagS3lossTo3Alt(param.binsInCoagulation,param.rho,param.T,param.visc,cmd,ln2s)*y(4)^2*1e6+ ...
                    coagTransferSFrom0To3;

                dy(6) = dy(6) + coagM3gainFrom0(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc,cmd,ln2s)*y(1)*y(4)*1e6 + ...
                    coagM3gainFrom3Alt(param.binsInCoagulation,param.rho,param.T,param.visc,cmd,ln2s)*y(4)^2*1e6 + ...
                    coagM3lossTo0(param.binsInCoagulation,alpha0,d1_0,d2_0,param.rho,param.T,param.visc,cmd,ln2s)*y(1)*y(4)*1e6 + ...
                    coagM3lossTo3Alt(param.binsInCoagulation,param.rho,param.T,param.visc,cmd,ln2s)*y(4)^2*1e6+ ...
                    coagTransferMFrom0To3;

            elseif param.modes == 3
                dy(4) = dy(4) + coagN3Alt(param.binsInCoagulation,cmd,ln2s,param.rho,param.T,param.visc)*y(4)^2*1e6;
                dy(5) = dy(5) + coagS3Alt(param.binsInCoagulation,cmd,ln2s,param.rho,param.T,param.visc)*y(4)^2*1e6;
            end
        end
    end
    
    
    % coagulational losses to bg
    if param.coagSink == 1
        dy(1) = dy(1) - param.coagSinkFactor*d1_0^param.coagSinkExponent*lossN0(alpha0,dd0,param.coagSinkExponent)*y(1)*param.coagSinkN*1e6;
        dy(2) = dy(2) - param.coagSinkFactor*d1_0^param.coagSinkExponent*lossS0(alpha0,dd0,param.coagSinkExponent)*y(2)*param.coagSinkN*1e6;
        dy(3) = dy(3) - param.coagSinkFactor*d1_0^param.coagSinkExponent*lossM0(alpha0,dd0,param.coagSinkExponent)*y(3)*param.coagSinkN*1e6;
        dy(4) = dy(4) - param.coagSinkFactor*cmd^param.coagSinkExponent*exp(param.coagSinkExponent^2/2*ln2s)*y(4)*param.coagSinkN*1e6;
        dy(5) = dy(5) - param.coagSinkFactor*cmd^param.coagSinkExponent*exp((param.coagSinkExponent^2/2+param.coagSinkExponent*2)*ln2s)*y(5)*param.coagSinkN*1e6;
        dy(6) = dy(6) - param.coagSinkFactor*cmd^param.coagSinkExponent*exp((param.coagSinkExponent^2/2+param.coagSinkExponent*3)*ln2s)*y(6)*param.coagSinkN*1e6;
    end

end
