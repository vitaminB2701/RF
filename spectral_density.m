function J = spectral_density(x,S0)
%SPECTRAL_DENSITY Calculate spectral density
%   Synthax
%       J = spectral_density(S)
%   
%   Description
%       Calculate the normalized spectral density function J0(x).
%       Müh et al 2010 JPC
%       Raszewski & Renger 2008 JACS (eg. 40)
%
%   Input
%     x - angular frequency [rad/s]
%     S0 - Huang-Rhys factor (default 1)
    
    if ~exist('S0','var')
        S0 = 1;
    end

    c = 2.997e-2; % speed of light in cm/ps

    s1 = 0.8; hw1 = 0.56; w1 = 2*pi*c*hw1; 
    s2 = 0.5; hw2 = 1.94; w2 = 2*pi*c*hw2;

    negx = x < 0;
    
    Jterm = @(s,w) s./(5040*2*w^4).*abs(x).^3.*exp(-sqrt(abs(x)./w));
    J = S0/(s1+s2) .* (Jterm(s1,w1) + Jterm(s2,w2));        
%     J(negx) = -J(negx);
  J(negx) = 0;
end

