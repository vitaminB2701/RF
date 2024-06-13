function K = genfoerster(X,E,U,V,Da,Di,ig,T)
%GENFOERSTER Calculate rates using generalized Förster theory
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
c0 = 2.997e10; % speed of light [cm/s]
kB = 0.695; % boltzmann constant [cm^-1/K]
kT = kB*T;  % boltzmann energy [cm^-1]

N = size(U,1);
K = zeros(N);

% Loop over excitonic clusters (domains)
Ng = size(ig,2);
for a = 1:Ng
    ia = ig(:,a);
    for b = 1:Ng
        if b~=a
            ib = ig(:,b);
            for m = find(ia)'
                for n = find(ib)'
                    
                    % Calculate coupling energy Vmn
                    % Raszewski & Renger 2008 eq. 8
                    if E(m) >= E(n) % Only for downhill transfer
%                         Vmn = U(ia,m)'*V(ia,ib)*U(ib,n); %This formula
%                         can remove the for loop for ma&nb, but somehow
%                         takes longer time
                        Vmn = 0;
                        for ma = 1:N
                            for nb = 1:N
                                Vmn = Vmn + U(ma,m)*U(nb,n)*V(ma,nb);
                            end
                        end
                        
                        % Calculate rate constant
                        % Raszewski & Renger 2008 eq. 21
                        % Pullerits et al 1997 JPC
                        K(n,m) = 1.183*Vmn^2*trapz(X,Di(:,m).*Da(:,n));
%                         K(n,m) = 1.183*Vmn^2*(trapz(X,Di(:,m).*Da(:,n))+trapz(X,Di(:,n).*Da(:,m)));
                        
                        % Calculate uphill rate using detailed balance
                        K(m,n) = exp(-(E(m)-E(n))/kT)*K(n,m);
%                         % Scale to detailed balance
%                         K(m,n) = 1/(1+exp(-(E(n)-E(m))/kT))*K(n,m);
%                         K(n,m) = 1/(1+exp(-(E(m)-E(n))/kT))*K(n,m);

                    end
                end
            end
        end
    end
end
% % Convert K to [ps^-1]
% K = K*1e-12;

end

