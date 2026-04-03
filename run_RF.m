% Run Redfield-Forster model
clear; %close all

%% Load input
input_RF;

%% Set up
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

%% Plot results
plot_RF;
