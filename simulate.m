% Runs a simulation

disp('Simulating...')

switch(p.model(1:2))
    case 'FS'
        out = simulationFS(p);
    otherwise
        alphaPrev=0;
        dPrev=0.5;
        out = simulationPLLN(p);
        clear alphaPrev dPrev
end

disp('Done')

if p.plotOutputAfterSim
    plotOutput(out)
end