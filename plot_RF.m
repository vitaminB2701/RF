%% Plot results

% Plot abs
figure; hold on;
plot(Par.X,RF.A./max(RF.A),Par.X,RF.F./max(RF.F)); 
xlabel('Wavenumber (cm-1)'); legend ('Abs','Fluo');
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
dat2D.Abs = rot90(permute(-RF.TA,[3 2 1]),2)+rot90(permute(-RF.TF,[3 2 1]),2);
%save 2D.mat dat2D
plot2Ds(dat2D,10,'unit','cm-1')