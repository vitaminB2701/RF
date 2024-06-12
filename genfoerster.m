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
                        V_AB = 0;
                        V_Anb = 0;
                        V_maB = 0;
                        k2 = 0; k3 = 0;
                        % Site loop
                        for nb = 1:N
                            for ma = 1:N
%                                 V_AB = V_AB + U(ma,A)*U(nb,B)*V(ma,nb)*vib.FC00sq;
                                V_Anb = V_Anb + U(ma,A)*V(ma,nb)*vib.FC00*vib.FC01;
                            end
                            k2 = k2 + 1.183*U(nb,B)^2*V_Anb^2*trapz(X,Di(:,A).*Dma(:,nb));
                        end
                        for ma = 1:N
                            for nb = 1:N
                                V_maB = V_maB + U(nb,B)*V(nb,ma)*vib.FC00*vib.FC01;
                            end
                            k3 = k3 + 1.183*U(ma,A)^2*V_maB^2*trapz(X,Dmi(:,ma).*Da(:,B));
                        end
                        
                        % Calculate rate constant
                        % Renger et al J Plan Physiol 2011
                        % Pullerits et al 1997 JPC
                        K(B,A) = 1.183*V_AB^2*trapz(X,Di(:,A).*Da(:,B)) + k2 + k3;
                        
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

