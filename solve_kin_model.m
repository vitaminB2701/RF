function [P, TA, TF, tau, a] = solve_kin_model(X,Xexc, t, K, Ae, Fe)
% Solve kinetic model and calculate species populations and fluorescence 

% Initial population
N = size(K,1);

%% Solve system
% Transfer matrix
K = K - diag(diag(K)) - diag(sum(K));

% Eigenvectors and eigenvalues
[U,L] = eig(K);
U1 = U';
tau = 1./diag(-L);
exptau = exp(-t./tau);

% Loop over excitations
P = zeros(N,numel(t),numel(Xexc));
TA = zeros(numel(t),numel(X),numel(Xexc));
TF = zeros(numel(t),numel(X),numel(Xexc));
for k = 1:numel(Xexc)
    iex = ind(X, Xexc(k)); % Excitation
    B = Ae(iex,:)'; % B is p(0)
%     B = B/sum(B);

    % Weighted eigenvectors
    a = (U\B).*U1;

    % Population kinetics
    p = real(a'*exptau);
    P(:,:,k) = p;

    % Absorption and fluorescence kinetics
    if exist('Ae','var')
        ta = Ae*p; TA(:,:,k) = ta';    
    end
    if exist('Fe','var')
        tf = Fe*p; TF(:,:,k) = tf';
    end
end