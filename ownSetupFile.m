p.condensationalTransfer = 1;
p.T = 280;
p.rho = 1400;
p.JMatrix = [0:100:18000; 0:0.01:1.8];
p.model = 'FS40';
p.initialMomentVec = zeros(1,n_sec(p.model));