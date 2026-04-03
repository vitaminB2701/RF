%% Input Parameters
Par.X = 14000:20:16000; % Wavenumber scale [cm-1]
Par.T = 300;         % Temperature [K]
Par.Vc = 20;        % Coupling cutoff (for clustering)
Par.Ec = 300;       % Energy difference cutoff (for clustering)
Par.Rc = 0.5;       % Protein correlation radius [nm]
Par.t = logspace(-2,3,150);  % Population time (for kinetics) [ps]
Par.exc = 1400:20:16000;    % Excitation frequency (for kinetics) [cm-1]
Par.BlockSize = 8;  % Number of parallel runs in a block, for static disorder
Par.Niter = 125; % Number of blocks to run (total=BlockSize*Niter)
Par.taudeph = 0.150; % Pure dephasing time (ps)
Par.energyfile = fullfile('Energy','C1S1.txt'); % File containing ID, site E,...
Par.pdbfile = {fullfile('pdb','5xnm.pdb'),'ABCDGYNS'}; % File containing pdb name and specify chains
% Mask to manipulate K matrix
Kmask = ones(length(Par.C.chain));
    % Use this section to manipulate Kmask, comment if unused
    %--------------------------------------------------------
% Block all connections between 2 groups of chains
% Kmask(~~sum(Par.C.chain=='4',2),~~sum(Par.C.chain=='123',2)) = 0.1; 
% Modify all connections from a group of chains
% Kmask(~~sum(Par.C.chain=='4R',2),~sum(Par.C.chain=='4R',2)) = 0;
    %--------------------------------------------------------
Par.Kmask = Kmask.*Kmask';

% Output file
fileout = "RF_out"+'_'+datestr(now,'yyyymmdd_HHMMss')+'.mat';