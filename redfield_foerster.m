function RF = redfield_foerster(atom,Epar,Par,vib,fileout)
% Redfield-F�rster model calculation
%
% Input
%   atom - Mg, Nb, Nd atom coordinates for each Chl
%   Epar - struct with site energies (E), Huang-Rhys factors (S) and
%          inhomogeneous broadening (cinh)
%   Par  - parameter struct
%
% Output
%   RF - struct with fields:
%   V  - coupling energies
%   G  - exciton clusters
%   ig - exciton cluster index (sites in rows, clusters in columns)
%   E  - exciton energies
%   U  - eigenvectors
%   K  - transfer rates
%   A  - absorption spectrum
%   F  - emission spectrum
%   Ae - exciton absorption spectra
%   Fe - exciton emission spectra
%   P  - exciton population kinetics

%% Parameters
T = Par.T;        % Temperature
X = Par.X;        % Wavenumber scale
Rc = Par.Rc;      % Protein correlation radius [nm]
Xexc = Par.exc;   % excitation frequency
t  = Par.t;       % population time

% Disorder averaging iterations
BlockSize = Par.BlockSize;  % Parallel runs in a block
Niter = Par.Niter;          % Number of blocks to average

% Exciton parameters
E0 = Epar.E;    % Site energies [cm^-1]
N = length(E0); % Number of Chls
S0 = Epar.S;    % Huang-Rhys factors
cinh = Epar.cinh; % inhomogeneous width [cm^-1]
kd = Epar.kd;   % diagonal decay rates
kd = diag(kd);
taudeph = Par.taudeph;

% Save file
if ~exist('fileout','var')
    fileout = 'RF_out';
end

% Constants
c0 = 2.997e-2; % speed of light [cm/ps]
kB = 0.695; % boltzmann constant [cm^-1/K]

% Convert wavenumbers to angular frequency [rad/s]
ang_freq = @(x) 2*pi*c0*x;

% Frequency axis
W = ang_freq(X);  % [rad/s]
numX = numel(X);
numExc = numel(Xexc); % number of excitation frequencies

%% Create clustered Hamiltonian
C = Par.C; % Coupling info
V = C.V;        % Coupling energies [cm^-1]
RF.V = V;
R = C.R;        % Center-to-center distances [nm]
D = C.D;        % Dipole strength [Debye^2]
Dvec = C.Dvec;  % Dipole vectors 

% Divide into clusters by coupling with a cutoff Vc
[G,ig] = cluster_by_coupling(V,Par.Vc,E0,Par.Ec);
RF.G = G; RF.ig = ig;

% Display clusters
if isfield(C,'chain') && isfield(C,'resname')
    RF.chain = C.chain;
    clusters = unique(G);
    for c = clusters'
        fprintf('\n%d: %d',c)
        ci = find(G==c);
        for i = ci'
            fprintf('%s ',[C.chain(i) '.' C.resname{i} num2str(C.molid(i))]);
        end
    end
    fprintf('\n');
end

%% Exciton-vibrational line broadening term
x = 1:1:600;    % [cm^-1]
w = ang_freq(x);% [rad/ps]

% Spectral density
J = spectral_density(w);

% Reduced reorganization energy (Er = lambda/S) [rad/ps];
Er0 = trapz(w,w.*J); 

% Number of vibrational quanta (Bose-Einstein distribution)
nvib = 1./(exp(x/(kB*T))-1); % (eq. S5)

% Time domain
t1 = (0:10:8000) * 1e-3;  % [ps]
nt = numel(t1);
it0 = ind(t1,0);
% Correlation function G(t)
Gt = zeros(1,nt);
for ti = 1:nt
    % (eq. S3)
    Gtw = (1+nvib).*J.*exp(-1i*w*t1(ti)) + nvib.*J.*exp(1i*w*t1(ti));
    Gt(ti) = trapz(w,Gtw);
end
RF.Gt = Gt;
% Recalculation with negative frequencies
x = -600:1:600;
w = ang_freq(x);

nvib = 1./(exp(x/(kB*T))-1); % (eq. S5)
nvibn = 1./(exp(-x/(kB*T))-1); % (eq. S5)

J = spectral_density(w);
Jn = spectral_density(-w);

%% Monte-Carlo sampling of disorder
rng;

fun_redfield_lineshapes = @redfield_lineshapes;
fun_linear_spectra = @linear_spectra;

fprintf('\nRunning...\n');
tic

% Initialize outputs
RF.E = zeros(N,1);  % Averaged exciton energies
RF.U = zeros(N,N);  % Averaged eigenvectors
RF.K = zeros(N,N);  % Averaged rate constants
RF.A = zeros(numX,1); % Averaged Absorption
RF.F = zeros(numX,1); % Averaged Fluorescence
RF.Ae = zeros(numX,N);
RF.Fe = zeros(numX,N);
RF.P  = zeros(N,numel(Par.t),numExc); % Averaged population
RF.TA1 = zeros(numel(Par.t),numX,numExc); % Time-resolved absorption
RF.TA2 = RF.TA1;
RF.TF1 = zeros(numel(Par.t),numX,numExc); % Time-resolved fluorescence
RF.TF2 = RF.TF1;

for iter = 1:Niter
    Ed = zeros(N,BlockSize);    % Averaged exciton energies
    Ud = zeros(N,N,BlockSize);  % Averaged eigenvectors
    Kd = zeros(N,N,BlockSize);  % Averaged rate constants
    Ad = zeros(numX,BlockSize); % Averaged Absorption
    Fd = zeros(numX,BlockSize); % Averaged Fluorescence
    Aed = zeros(numX,N,BlockSize); % Exciton absorption
    Fed = zeros(numX,N,BlockSize); % Exciton fluorescence
    Pd = zeros(N,numel(Par.t),numExc,BlockSize); % Averaged population
    TA1d = zeros(numel(Par.t),numX,numExc,BlockSize); % Time-resolved absorption
    TA2d = TA1d;
    TF1d = zeros(numel(Par.t),numX,numExc,BlockSize); % Time-resolved fluorescence
    TF2d = TF1d;    
    
    for bl = 1:BlockSize
        % Random site energies
        Em = E0 - randn(N,1).*cinh;
        
        % Redfield and exciton lineshapes
        [E,U,Kr,Da,Di,Dma,Dmi] = feval(fun_redfield_lineshapes,Em,vib);
        
        % F�rster rate constants
        Kf = genfoerster(X,E,U,V,Da,Di,ig,T,vib,Dma,Dmi);
        
        % Save
        Ed(:,bl) = E;
        Ud(:,:,bl) = U;
        K = Kr + Kf + kd;
        Kd(:,:,bl) = K;
        
        % Calculate linear spectra
        [Ad(:,bl),Ae,Fe,mux] = feval(fun_linear_spectra,U,Da,Di,Dma,Dmi,vib);
        Aed(:,:,bl) = Ae;
        Fed(:,:,bl) = Fe;
        
        % Calculate polarization scaling factor
        mue = mux./sqrt(sum(mux.^2,2)); % Exciton dipole moment (unit vector)
        leg2 = (3*(mue*mue').^2-1)/2; % Legendre poly2 for each pair of excitons
        pol_par = (4*leg2+5)/45; % Factor for parallel pumpprobe
        pol_per = (5-2*leg2)/45; % Factor for perpendicular pumpprobe
        
        % Kinetics
        [P1,P2,TA1d(:,:,:,bl),TA2d(:,:,:,bl),TF1d(:,:,:,bl),TF2d(:,:,:,bl)]...
            = solve_kin_model(X,Xexc,t,K,Ae,Da,Di,pol_par,pol_per,mux);        
        Pd(:,:,:,bl) = P1+2*P2;
        
        % Steady-state emission
        p_i = trapz(t,mean(P1+2*P2,3)');
        Fd(:,bl) = sum(Fe.*p_i,2);
    end
    
    % Sum
    ri = (1-1/iter);
    RF.E = RF.E*ri + mean(Ed,2)/iter;
    RF.U = RF.U + mean(Ud,3);
    RF.U = RF.U ./ sqrt(sum(RF.U.^2,1));
    RF.K = RF.K*ri + mean(Kd,3)/iter;
    RF.A = RF.A*ri + mean(Ad,2)/iter;
    RF.F = RF.F*ri + mean(Fd,2)/iter;
    RF.Ae = RF.Ae*ri + mean(Aed,3)/iter;
    RF.Fe = RF.Fe*ri + mean(Fed,3)/iter;
    RF.P  = RF.P*ri + mean(Pd,4)/iter;
    RF.TA1 = RF.TA1*ri + mean(TA1d,4)/iter;
    RF.TA2 = RF.TA2*ri + mean(TA2d,4)/iter;
    RF.TF1 = RF.TF1*ri + mean(TF1d,4)/iter;
    RF.TF2 = RF.TF2*ri + mean(TF2d,4)/iter;
    
    RF.Efull(:,(iter-1)*BlockSize+1:iter*BlockSize) = Ed;
    RF.Ufull(:,:,(iter-1)*BlockSize+1:iter*BlockSize) = Ud;
    RF.Kfull(:,:,(iter-1)*BlockSize+1:iter*BlockSize) = Kd;
    
    % Save result by progress
    elapsed_time = toc;
    save(fileout,'atom','Epar','Par','X','RF','iter','C','elapsed_time')
    
    % Display progress
    if iter==1, fprintf('\nIterations:'); end
    if rem(iter,20)==1, fprintf('\n'); end
    fprintf(' %d',iter)
    
end
fprintf('\nDone!\n');


%% Subroutines
    % Calculate spectra and dynamics for one realization
    function [E,U,Kr,Da,Di,Dma,Dmi] = redfield_lineshapes(Em,vib)
        % Diagonalize partitioned Hamiltonian
        [E,U] = diag_hamiltonian(Em,vib);
        
        Kr = zeros(N); % Redfield rates
        Yr = eye(N); % Redfield Y
        
        Da = zeros(length(X),N); % absorption lineshape
        Di = zeros(length(X),N); % emission lineshape
        [Dma,Dmi] = intraviblineshape(Em,vib); % localized intravibronic lineshape
        for k = 1:size(ig,2)
            ix = ig(:,k); % cluster indices
            
            % Exciton Huang-Rhys factors
            Uc = U(ix,ix);
            Sc = Uc'.^2 * S0(ix);
            
            % Calculate Redfield rates
            if numel(find(ix))>1
                [Kr(ix,ix), Yr(ix,ix)] = redfield(E(ix),Uc,R(ix,ix),Rc,Sc,T);
            end
            
            % Calculate lineshapes
            [Da(:,ix),Di(:,ix)] = lineshape(E(ix),Sc,Kr(ix,ix),Yr(ix,ix));
        end
    end

    % Diagonalize Hamiltonian by clusters
    function [E,U] = diag_hamiltonian(Em,vib)
        n = numel(Em);
        H = diag(Em) + V*vib.FC00sq;
        E = zeros(n,1);
        U = zeros(n);
        for g = 1:size(ig,2)        
            ix = ig(:,g);
            Hp = H(ix,ix);
            [Up,Ep] = eig(Hp);
            E(ix) = diag(Ep);
            U(ix,ix) = Up;
        end
    end

    %Calculate absorption and emission exciton lineshapes
    function [Dma,Dmi] = intraviblineshape(Em,vib)
        Dma = zeros(numX,numel(Em));
        Dmi = zeros(numX,numel(Em));
        % Transition dipole moments
        for nChl = 1:numel(Em)
            intraviba = zeros(1,numel(t1));
            intravibi = zeros(1,numel(t1));
            for nvibmode = 1:numel(vib.vibmode) %intravib modes
                wma = ang_freq(Em(nChl)+vib.vibmode(nvibmode))-Er0*S0(nChl);
                wmi = ang_freq(Em(nChl)-vib.vibmode(nvibmode))-Er0*S0(nChl);
                intraviba = intraviba + vib.FCi01sq(nvibmode)*exp(1i.*(W'-wma)*t1).*exp(S0(nChl)*(Gt-Gt(1))-abs(t1)/taudeph);
                intravibi = intravibi + vib.FCi01sq(nvibmode)*exp(-1i.*(W'-wmi)*t1).*exp(S0(nChl)*(Gt-Gt(1))-abs(t1)/taudeph);
            end
            Dma(:,nChl) = real(trapz(t1,intraviba,2))/(2*pi);
            Dmi(:,nChl) = real(trapz(t1,intravibi,2))/(2*pi);
        end
        % Normalize
%         Dma = Dma./trapz(X,Dma,1);
%         Dmi = Dmi./trapz(X,Dmi,1);
    end

    function [Da,Di] = lineshape(E,S,K,Y)
        % Raszewski & Renger 2008
        n = numel(E);
        
        %% Calculate exciton lineshapes
        Da = zeros(numX,n);
        Di = zeros(numX,n);
        
         % Loop over excitons
        for k = 1:n
            
            % Exciton relaxation lifetime            
%             tau = 2*1e-12./sum(K(:,k)); % [s] (eq. S6)
%             if tau > taudeph
%                 tau = taudeph;
%             end
            tau = 1/(0.5*sum(K(:,k)) + 1/taudeph);
            
            % Total broadening term
            eVib = exp(S(k).*Y(k,k)*(Gt-Gt(it0)) - abs(t1)/tau); 
            
            % Calculate renormalized 0-0 transition frequency
            
            % Muh et al 2010 JPC, eq. S14
            % Raszewski & Renger 2008, eq. S7
            ck = 0;
            for m = 1:n
                if m~=k                    
                    % Frequency difference between states
                    wkm = ang_freq(E(k)-E(m));
                    
                    % Correlation function
                    Ck = pi.*w.^2.*((1+nvib).*S(k).*J+nvibn.*S(k).*Jn)./(wkm-w);
                    
                    % "clean up" the singulatities to calculate
                    % "principal value" of the integral
                    good = ~(isnan(Ck) | isinf(Ck));
                    ckm = trapz(w(good),Ck(good))/pi;                    
                    
                    ck = ck + Y(k,m)*ckm;
                end
            end
            
            % Renormalized frequency
            wk = ang_freq(E(k)) - Y(k,k)*Er0*S(k) + ck;
            
            % Exciton lineshape
            a = exp(1i.*(W' - wk)*t1).*eVib;
            f = exp(-1i.*(W' - wk)*t1).*eVib;
            Da(:,k) = real(trapz(t1,a,2))/(2*pi);
            Di(:,k) = real(trapz(t1,f,2))/(2*pi);
                       
            % Normalize the lineshapes on a wavenumber scale
%             Da(:,k) = Da(:,k) ./ trapz(X,Da(:,k));
%             Di(:,k) = Di(:,k) ./ trapz(X,Di(:,k));
        end
    end

%   Calculate linear spectra
    function [A,Ae,Fe,mux] = linear_spectra(U,Da,Di,Dma,Dmi,vib)
        
        % Transition dipole moments
        n = 1.4;          % refractive index
        mu = sqrt(D*n).*Dvec;   % monomeric dipole moments
        mux = U'*mu*vib.FC00;       % excitonic dipole moments
        mux2 = sum(mux'.^2); % Square exciton dipole strengths
        Dma_x = (n*D'.*Dma)*(U.^2);
        Dmi_x = (n*D'.*Dmi)*(U.^2);
%         % Absorption
        Ae = Da.*mux2 + Dma_x;  % exciton absorption spectra
        A = sum(Ae,2); % absorption spectrum
        
        % Fluorescence
        Fe = Di.*mux2 + Dmi_x; % exciton emission spectra
        
        % kT = kB*T;  % boltzmann energy [cm^-1]
        % fB = exp(-E(:)/kT); fB = fB'/sum(fB);
        % F = sum(fB.*Fe,2);     % emission spectrum
    end
end