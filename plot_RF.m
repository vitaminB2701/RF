%% Plot results

% Plot abs
figure; hold on;
plot(1e7./Par.X,RF.A./max(RF.A),1e7./Par.X,RF.F./max(RF.F)); 
xlabel('Wavelength (nm)'); legend ('Abs','Fluo');
box on;

%% Coupling strengths and EET
if size(Epar,1)>2
    plot_fret_struct(atom,RF.V,RF.G,[1 100],num2cellstr(Epar.Chl),true)
    plot_fret_struct(atom,RF.K,RF.G,[0.1 5],num2cellstr(Epar.Chl),true)
end

%% Population
figure; formatfig('default');
exc_range = unit_conv([600 700],'nm to cm-1');
exci = ind(Par.exc,exc_range(2)):ind(Par.exc,exc_range(1));
for ch = unique(RF.chain)'
    plot(Par.t,sum(sum(RF.P(RF.chain==ch,:,exci),3)),'linewidth',1.2,'displayname',ch);
    hold on;
end
legend;
set(gca,'xscale','log','fontsize',14);
xlabel('Time (ps)'); ylabel('Population');
xlim([10^-2,10^3]);
box on;

%% Plot 2D
dat2D.X = flip(1e7./Par.X);
dat2D.Y = flip(1e7./Par.exc)';
dat2D.T = Par.t;
dat2D.Abs = rot90(permute(-RF.TA,[3 2 1]),2)+rot90(permute(-RF.TF,[3 2 1]),2);
%save 2D.mat dat2D
plot2Ds(dat2D,0)