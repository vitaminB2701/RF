% Redfield-Förster model
%
% Input
% 1) PDB file containing coordinates of Mg, Nb and Nd atoms for each Chl
% 2) Excel table with site energies

clear; %close all

%% Set-up
% Parameters
Par.X = 14000:20:16000; % Wavenumber scale
Par.T = 80;         % Temperature [K]
Par.Vc = 20;        % Coupling cutoff (for clustering)
Par.Ec = 300;       % Energy difference cutoff (for clustering)
Par.Rc = 0.5;       % Protein correlation radius [nm]
Par.t = logspace(-2,3,150);  % Population time (for kinetics) [ps]
Par.exc = 14400:20:16000;    % Excitation frequency (for kinetics)
Par.BlockSize = 4;  % Number of parallel runs in a block, for static disorder
Par.Niter = 25; % Number of blocks to run (total=BlockSize*Niter)
Par.taudeph = 0.3; % Pure dephasing time (ps)

% Output file
fileout = "RF_out"+'_'+datestr(now,'yyyymmdd_HHMMss')+'.mat';

% Import exciton parameter table
Epar = readtable('Energy.xlsx');

% Structure file
atom = import_pdb('5xnm.pdb','G'); % Load structure

% Pick a few chromophores
iPick = [1 2];
Epar = Epar(iPick,:);
[chList,chListIndex] = uniquearray([atom.molid]);
chList(chList>614) = [];
atomList = [];
for i = 1:length(iPick)
    atomList = [atomList chListIndex(iPick(i)):chListIndex(iPick(i)+1)-1];
end
atom = atom(atomList);

% Coupling information
Par.C = calc_coupling(atom);

%% Run calculation
if isempty(gcp('nocreate'))==1
    parpool('local', Par.BlockSize);
end
RF = redfield_foerster(atom,Epar,Par,fileout);