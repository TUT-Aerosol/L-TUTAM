p.condensationalTransfer = 1;
p.T = 280;
p.rho = 1400;
p.JMatrix = [0:100:3600; 0:0.01:.36];
p.model = 'FS40';
p.initialMomentVec = zeros(1,n_sec(p.model));
p.losses=0;
p.coagSink=0;