function [P1, P2, TA1, TA2, TF1, TF2, tau, a] = solve_kin_model(X,Xexc, t, K, Ae, Da, Di, pol_par, pol_per, mux)
% Solve kinetic model and calculate species populations and fluorescence 

% Initial population
N = size(K,1);

% Square exciton dipole strengths
mux2 = sum(mux'.^2);

%% Solve system
% Transfer matrix
K = K - diag(diag(K)) - diag(sum(K));

% Eigenvectors and eigenvalues
[U,L] = eig(K);
U1 = U';
tau = 1./diag(-L);
exptau = exp(-t./tau);


% Initialize spectra matrices (1: parallel, 2: perpendicular)
P1 = zeros(N,numel(t),numel(Xexc));
P2 = P1;
TA1 = zeros(numel(t),numel(X),numel(Xexc));
TA2 = TA1;
TF1 = zeros(numel(t),numel(X),numel(Xexc));
TF2 = TF1;

% Loop over excitations
for k = 1:numel(Xexc)
    iex = ind(X, Xexc(k)); % Excitation
    B = Ae(iex,:)'; % B is p(0)
    B = diag(B); % Turn B to diagonal matrix
    
    % Iterate with each exciton
    for iB = 1:length(B)
        Bx = B(:,iB);
        a = (U\Bx).*U1; % Weighted eigenvectors
        p = real(a'*exptau); % Population kinetics
        
        % Apply polarization effects
        p_par = pol_par(:,iB).*p;
        p_per = pol_per(:,iB).*p;
        
        % Combine with other excitons
        P1(:,:,k) = P1(:,:,k) + p_par;
        P2(:,:,k) = P2(:,:,k) + p_per;

    end

    % Absorption and fluorescence kinetics
    if exist('Da','var')
        ta1 = (mux2.*Da)*P1(:,:,k); TA1(:,:,k) = ta1';    
        ta2 = (mux2.*Da)*P2(:,:,k); TA2(:,:,k) = ta2';

    end
    if exist('Di','var')
        tf1 = (mux2.*Di)*P1(:,:,k); TF1(:,:,k) = tf1';
        tf2 = (mux2.*Di)*P2(:,:,k); TF2(:,:,k) = tf2';
    end
end