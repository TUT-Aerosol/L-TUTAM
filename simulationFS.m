function out =  simulationFS(p)
%SIMULATIONFS simulates the FS model using input parameter set 'p'    

    nSec = n_sec(p.model);
    Dp_edges = logspace(log10(p.dCluster),log10(p.highestDiameter),nSec+1);
    logEdges = log10(Dp_edges);
    logCenters = logEdges(1:end-1) + diff(logEdges)./2;
    Dp_centers = 10.^(logCenters); 
    p.DeltaDp = diff(Dp_edges);
    p.Dp_centers = Dp_centers;

    p.kk = zeros(nSec,nSec); % Preallocate coagulation coefficients
    if p.coag
        for i = 1:nSec
            p.kk(i,:) = 2*pi*coag_kernel_dp(Dp_centers(i),Dp_centers,p.T,p.rho)*1e6;
        end
    end  

    % Calculating coagulation sink factor
    p.coagSinkFactor = coag_kernel_dp(5e-9,p.coagSinkCMD,p.T,p.rho)*2*pi*5e-9^(-p.coagSinkExponent);

    % Total simulation time
    p.totalTime = p.timeVec(end)-p.timeVec(1);
    
    % If 'p.initialMomentVec' is scalar 0, make a vector
    if length(p.initialMomentVec)==1
        if p.initialMomentVec==0
            p.initialMomentVec = zeros(1,nSec);
        end
    end

    % Options to ODE solver
    options = odeset('RelTol',p.relativeTolerance,'nonnegative',1:nSec,'stats','off');
    
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
    out.t = t; 
    out.Y = Y;
    out.p = p;

    % Zeros
    out.N = zeros(length(t),1);
    out.M_2 = out.N;
    out.M_3 = out.N;
    out.GMD = out.N;
    out.GSD = out.N;

    % Update values
    for i=1:length(t)
        out.N(i) = sum(Y(i,:));
        out.M_2(i) = sum(Y(i,:).*p.Dp_centers.^2);
        out.M_3(i) = sum(Y(i,:).*p.Dp_centers.^3);
        lngmd = sum(Y(i,:).*log(p.Dp_centers))/sum(Y(i,:));
        ln2gsd = sum(Y(i,:).*(log(p.Dp_centers)-lngmd).^2)/sum(Y(i,:));
        out.GMD(i) = exp(lngmd);
        out.GSD(i) = exp(sqrt(ln2gsd));
    end

    % Removing a large interpolation table from the output struct
    out.p=rmfield(out.p,'alphaTaulukko');
    out.p=rmfield(out.p,'dTaulukko');

end



function dy = modelFunc(t,y,param)
%MODELFUNC is the internal function containing the FS model equations

    if param.plotWaitbarDuringSim
        waitbar(t/param.totalTime);
    end

    % Find J and GR from the matrices
    param.J=param.JMatrix(2,find(t<=param.JMatrix(1,:),1));
    param.GR=param.GRMatrix(2,find(t<=param.GRMatrix(1,:),1));
    
    % Plotting the distr.
    if param.plotDistrDuringSim
        figure(1)
        clf


        stairs(param.Dp_centers*1e9,y/(log10(param.Dp_centers(3))-log10(param.Dp_centers(2))),'r')
        set(gca,'xscale','log','yscale','log')
        ylim([max(param.JMatrix(2,:))/100 max(param.JMatrix(2,:))*10]*param.totalTime)
        xlim([1 100])
        title(strcat('t= ',num2str(t,'%1.3f'),' s'))
        xlabel('Dp (nm)')
        ylabel('dN/dlogDp (cm^{-3})')

        drawnow
    end
    

    dy = zeros(size(y));

    for i = 1:length(dy)
       % new particle formation
        if i==1
            dy(i) = dy(i) + param.J;
        end

        % condensation
        C = param.GR/param.DeltaDp(i);

        if i>1,
            dy(i) = dy(i) - y(i).*C + y(i-1).*C;
        else
            dy(i) = dy(i) - y(i).*C;
        end

        % depositional losses
        if param.losses == 1
            dy(i) = dy(i) - param.lossesCoeff/(param.Dp_centers(i))*y(i);
        end


        % coagulation
        if param.coag
                cM = coagulationMatrix(param.Dp_centers,i);
                
                for j = 1:i
                    NNkk = y(i)*y(j)*param.kk(i,j);
                    
                    if i == j
                        dy(j) = dy(j)-NNkk; % loss      
                        dy = dy + cM(j,:)'.*0.5.*NNkk; % gain
                    else
                        dy(i) = dy(i)-NNkk; %loss
                        dy(j) = dy(j)-NNkk; %loss           
                        dy = dy + cM(j,:)'.*NNkk; % gain
                    end
                end
        end



    % coagulational losses
    if param.coagSink
        dy(i) = dy(i) - param.coagSinkFactor*param.coagSinkN*y(i)*param.Dp_centers(i)^param.coagSinkExponent*1e6;    
    end


    end



end
