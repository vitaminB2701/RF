% Redfield-Forster model
%
% Input
% 1) PDB file containing coordinates of Mg, Nb and Nd atoms for each Chl
% 2) Excel table with site energies

% Static disorder is simulated by repeating the calculation with Gaussian
% randomized site energies. The sampling is divided into blocks, where
% after each block is done calculating, the result is averaged with
% previous blocks and saved. The process is iterated so that the total
% randomizations samples is Nsample=BlockSize*Niteration. 
% By this way, the process of sampling over static disorder can be
% constantly recorded throughout the calculation. The program can be 
% terminated at anytime during the sampling process without losing the 
% already simulated results.
% However, a portion of calculation time will be sacrificed to write data.
% If desired, the number of iteration (Niter) can be reduced while
% increasing the BlockSize to reduced the writing time.
% Generally, a parallel run does not improve the performance of
% this code.

clear; %close all

%% Set-up
% Parameters
Par.X = 14280:20:16700; % Wavenumber scale [cm-1]
Par.T = 300;         % Temperature [K]
Par.Vc = 20;        % Coupling cutoff (for clustering)
Par.Ec = 300;       % Energy difference cutoff (for clustering)
Par.Rc = 0.5;       % Protein correlation radius [nm]
Par.t = logspace(-2,3,150);  % Population time (for kinetics) [ps]
Par.exc = 14280:20:16700;    % Excitation frequency (for kinetics) [cm-1]
Par.BlockSize = 6;  % Number of parallel runs in a block, for static disorder
Par.Niter = 100; % Number of blocks to run (total=BlockSize*Niter)
Par.taudeph = 0.150; % Pure dephasing time (ps)
Par.energyfile = fullfile('Energy','LHCIIpent_cohen2009_250shifted.txt'); % File containing ID, site E,...
Par.pdbfile = {fullfile('pdb','5xnm_123R_CLA616removed.pdb'),'1234R'}; % File containing pdb name and specify chains

% Output file
fileout = "RF_out"+'_'+datestr(now,'yyyymmdd_HHMMss')+'.mat';

% Import exciton parameter table
Epar = readtable(Par.energyfile);

% Structure file
atom = import_pdb(Par.pdbfile{1},Par.pdbfile{2}); % Load structure

% Coupling information
Par.C = calc_coupling(atom);

% Check consistency between structure and energy file
if ~isequal(Par.C.molid,Epar.ID)
    error("Molecule IDs in structure and energy files do not match");
end

%% Run calculation
RF = redfield_foerster(atom,Epar,Par,fileout);
