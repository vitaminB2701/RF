function [K,Y] = redfield(E,c,R,Rc,S,T)
%REDFIELD Calculate Redfield exciton relaxation rate constant
%  Raszewski & Renger 2008
%   
% Synthax
%    [K,Y] = redfield(E,U,S,T,R,Rc)
% 
% Description
%    Calculates the matrix of relaxation rate constants K [ps^-1]
%    and the matrix of amplitudes Y.
%
% Input
%   E - exciton energies (eigenstates) [cm^-1]
%   U - eigenvectors
%   R - center-to-center distances between sites [nm]
%   Rc - correlation radius of protein vibrations [nm]
%   S - Huang-Rhys factors
%   T - temperature [K]


% Constants
c0 = 2.997e-2; % speed of light [cm/ps]
kB = 0.695; % boltzmann constant [cm^-1/K]
kT = kB*T;  % boltzmann energy [cm^-1]

N = numel(E); % Number of states
if isscalar(S)
    S = repmat(S,N,1);
end

% Convert wavenumbers to angular frequency [rad/s]
ang_freq = @(x) 2*pi*c0*x; 

%% Calculate Y matrix
% Loop over exciton states k and l
Y = zeros(N);
for k = 1:N
    for l = k:N        
        y = 0;
        % Loop over sites m and n
        for m = 1:N
            for n = 1:N
                y = y + c(m,k)*c(m,l)*c(n,k)*c(n,l)*exp(-R(m,n)/Rc);
            end
        end
        Y(k,l) = y;  
        Y(l,k) = y;  
    end
end

%% Calculate rates
K = zeros(N,N);
for k = 1:N
    for l = 1:N
        if k~=l && E(k) >= E(l)
            % Transition frequency
            dE = E(k)-E(l);
            dw = ang_freq(dE);
            
            % Spectral density at that frequence            
            Jk = spectral_density(dw,S(k));
            Jl = spectral_density(-dw,S(l));

            % Number of vibrational quanta
            nk = 1./(exp(dE/kT)-1);
            nl = 1./(exp(-dE/kT)-1);
            
            % Rate constants
            K(l,k) = 2*pi*Y(k,l)*dw^2*(Jk*(nk+1) + Jl*nl);
            K(k,l) = 2*pi*Y(l,k)*dw^2*(Jl*(nl+1) + Jk*nk);
        end
    end
end

end

