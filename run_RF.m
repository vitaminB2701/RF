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
