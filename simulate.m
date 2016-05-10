% Runs a simulation

disp('Simulating...')

switch(p.model(1:2))
    case 'FS'
        out = simulationFS(p);
    otherwise
        out = simulationPLLN(p);
end

disp('Done')

if p.plotOutputAfterSim
    plotOutput(out)
end