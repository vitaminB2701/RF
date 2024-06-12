%% Plot results

% Plot abs
figure; hold on;
plot(Par.X,RF.A./max(RF.A),Par.X,RF.F./max(RF.F)); 
xlabel('Wavenumber (cm-1)'); legend ('Abs','Fluo');
box on;

% Plot abs in wavelegnth
figure; hold on;
plot(1e7./Par.X,RF.A./max(RF.A),1e7./Par.X,RF.F./max(RF.F)); 
xlabel('Wavelength (nm)'); legend ('Abs','Fluo');
box on;

%% Coupling strengths and EET
if size(Epar,1)>2
    plot_fret_struct(atom,RF.V,RF.G,[1 100],num2cellstr(Epar.ID),true)
    title('Couling map');
    plot_fret_struct(atom,RF.K,RF.G,[0.1 5],num2cellstr(Epar.ID),true)
    title('EET map');
end

%% Plot 2D
dat2D.X = flip(Par.X);
dat2D.Y = flip(Par.exc)';
dat2D.T = Par.t;
dat2D.Abs1 = rot90(permute(-RF.TA1,[3 2 1]),2)+rot90(permute(-RF.TF1,[3 2 1]),2);
dat2D.Abs2 = rot90(permute(-RF.TA2,[3 2 1]),2)+rot90(permute(-RF.TF2,[3 2 1]),2);
dat2D.Abs = 1/3*dat2D.Abs1 + 2/3*dat2D.Abs2;
%save 2D.mat dat2D
plot2Ds(dat2D,1,'unit','cm-1')