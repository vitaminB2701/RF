% Redfield-Forster model
%
% Input
% 1) PDB file containing coordinates of Mg, Nb and Nd atoms for each Chl
% 2) Excel table with site energies

clear; %close all

%% Set-up
% Parameters
Par.X = 14000:20:16000; % Wavenumber scale [cm-1]
Par.T = 300;         % Temperature [K]
Par.Vc = 20;        % Coupling cutoff (for clustering)
Par.Ec = 300;       % Energy difference cutoff (for clustering)
Par.Rc = 0.5;       % Protein correlation radius [nm]
Par.t = logspace(-2,3,150);  % Population time (for kinetics) [ps]
Par.exc = 14000:20:16000;    % Excitation frequency (for kinetics) [cm-1]
Par.BlockSize = 6;  % Number of parallel runs in a block, for static disorder
Par.Niter = 100; % Number of blocks to run (total=BlockSize*Niter)
Par.taudeph = 0.150; % Pure dephasing time (ps)
Par.energyfile = fullfile('Energy','LHCIItri.txt'); % File containing ID, site E,...
Par.pdbfile = {fullfile('pdb','5xnm.pdb'),'GNY'}; % File containing pdb name and specify chains

% Output file
fileout = "RF_out"+'_'+datestr(now,'yyyymmdd_HHMMss')+'.mat';

% Import exciton parameter table
Epar = readtable(Par.energyfile);

% Structure file
atom = import_pdb(Par.pdbfile{1},Par.pdbfile{2}); % Load structure

% Coupling information
Par.C = calc_coupling(atom);

%% Run calculation
% Initiate parallel pool if no pool presents
%if isempty(gcp('nocreate'))==1
%    parpool('local', Par.BlockSize);
%end
RF = redfield_foerster(atom,Epar,Par,fileout);
