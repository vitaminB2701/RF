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
Par.T = 80;         % Temperature [K]
Par.Vc = 20;        % Coupling cutoff (for clustering)
Par.Ec = 300;       % Energy difference cutoff (for clustering)
Par.Rc = 0.5;       % Protein correlation radius [nm]
Par.t = logspace(-2,3,150);  % Population time (for kinetics) [ps]
Par.exc = 14280:20:16700;    % Excitation frequency (for kinetics) [cm-1]
Par.BlockSize = 6;  % Number of parallel runs in a block, for static disorder
Par.Niter = 20; % Number of blocks to run (total=BlockSize*Niter)
Par.taudeph = 0.300; % Pure dephasing time (ps)
Par.energyfile = fullfile('Energy','Chlab.txt'); % File containing ID, site E,...
Par.pdbfile = {fullfile('pdb','Chlab.pdb'),'Y'}; % File containing pdb name and specify chains

% Output file
fileout = "RF_out"+'_'+datestr(now,'yyyymmdd_HHMMss')+'.mat';

% Import exciton parameter table
Epar = readtable(Par.energyfile);

% Import intravibration parameters
vibpar = load('Vibpar.txt');
vib.vibmode = vibpar(:,1); % Vibrational modes [cm-1]                                       
vib.vibHR = vibpar(:,2);  % HuangRhys factors
scale_vib = nthroot(prod(exp(-vib.vibHR)),length(vibpar)); % To normalize FC00 to 1
% Frank-Condon factors
vib.FCi00sq = exp(-vib.vibHR)/scale_vib; 
vib.FCi01sq = exp(-vib.vibHR).*vib.vibHR/scale_vib;
vib.FCi00 = sqrt(vib.FCi00sq);
vib.FCi01 = sqrt(vib.FCi01sq);
vib.FC00sq = prod(vib.FCi00sq);
vib.FC00 = sqrt(vib.FC00sq);
vib.FC01sq = prod(vib.FCi01sq);
vib.FC01 = sqrt(vib.FC01sq);

% Structure file
atom = import_pdb(Par.pdbfile{1},Par.pdbfile{2}); % Load structure

% Coupling information
Par.C = calc_coupling(atom);

%% Run calculation
RF = redfield_foerster(atom,Epar,Par,fileout);
