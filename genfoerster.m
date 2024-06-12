function K = genfoerster(X,E,U,V,Da,Di,ig,T,vib,Dma,Dmi)
% GENFOERSTER Calculate rates using generalized Förster theory
%  Raszewski & Renger 2008, eq. 21
%
% Synthax
%   K = genforster(X,Da,Di,ig)
%
% Input
%   X  - wavenumber scale
%   U  - exciton eigenvectors
%   V  - site couplings
%   Da - exciton absorption lineshapes (see lineshape)
%   Di - exciton emission lineshapes
%   ig - exciton cluster map (see cluster_by_coupling)
%

% Constants
c0 = 2.998e10; % speed of light [cm/s]
kB = 0.695; % boltzmann constant [cm^-1/K]
kT = kB*T;  % boltzmann energy [cm^-1]

N = size(U,1);
K = zeros(N);

% Area-normalizing exciton lineshapes and intravib lineshapes
Da = Da./trapz(X,Da,1);
Di = Di./trapz(X,Di,1);
if sum(vib.vibHR)~=0
    Dma = Dma./trapz(X,Dma,1);
    Dmi = Dmi./trapz(X,Dmi,1);
end

% Loop over excitonic clusters (domains)
Ng = size(ig,2);
% Domain loop
for a = 1:Ng
    ia = ig(:,a);
    for b = 1:Ng
        if b~=a
            ib = ig(:,b);
            % Exciton loop
            for A = find(ia)'
                for B = find(ib)'                    
                    % Calculate couplings V_MN V_Mnb V_maN
                    % Renger et al J Plan Physiol 2011
                    if E(A) >= E(B) % Only for downhill transfer
                        % Exciton-exciton coupling
                        V_AB = U(:,A)'*V*U(:,B)*vib.FC00sq;
                        % Exciton-vibration coupling
                        V_Anb = (U(:,A)'*V(:,ib))'.*U(ib,B)*vib.FC00*vib.FC01;
                        V_maB = (U(:,B)'*V(:,ia))'.*U(ia,A)*vib.FC00*vib.FC01;
                        % Rate components
                        k1 = 1.183*trapz(X,Di(:,A).*Da(:,B))*V_AB^2;
                        k2 = 1.183*trapz(X,Di(:,A).*Dma(:,ib))*(V_Anb.^2);
                        k3 = 1.183*trapz(X,Dmi(:,ia).*Da(:,B))*(V_maB.^2);

                        % Calculate rate constant
                        % Renger et al J Plan Physiol 2011
                        % Pullerits et al 1997 JPC
                        K(B,A) = k1 + k2 + k3;
                        
                        % Calculate uphill rate using detailed balance
                        K(A,B) = exp(-(E(A)-E(B))/kT)*K(B,A);
                    end
                end
            end
        end
    end
end
% Note: the 1.183 factor is to convert trapz (intergral) in cm-1, couplings
% in cm-1 into 2pi/hbar in SI and yield K in 1/ps
end

