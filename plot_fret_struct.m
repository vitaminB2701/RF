function plot_fret_struct(atom,K,G,krange,names,plot3D)
%PLOT_FRET_STRUCT Overlay FRET connections on the molecular structure

if ~exist('plot3D','var')
    plot3D = false;
end

mg = atom(strcmp({atom.type},'MG'));
N = length(mg);
molid = [mg.molid];

if ~exist('krange','var')
    krange = [10 1000];
end
kmin = krange(1); kmax = krange(2);
thickness = 2;

if ~exist('names','var')
    names = num2str(molid(m));
end

%% Plot molecular structure
figure; hold on
cs = colorscheme('default');
% plotc = @(x,y,z) plot3(x,y,z,'o','markersize',8,...
%     'markerfacecolor',cs(get(gca,'colororderindex'),:));
% splitapply(plotc, [mg.x]', [mg.y]', [mg.z]', G)

%% Overlay rates
for m = 1:N
    for n = m+1:N
        k = K(m,n) + K(n,m);
        if k > kmin
            lw = ((log10(k)-log10(kmin))+0.2) .* thickness;
            if k < kmax
                clr = repmat((kmax-k)/kmax,1,3);
            else
                clr = [0 0 0];
            end
            if plot3D
                line([mg(m).x,mg(n).x],[mg(m).y,mg(n).y],[mg(m).z,mg(n).z],...
                'LineWidth',lw,'Color',clr)
            else
                line([mg(m).x,mg(n).x],[mg(m).z,mg(n).z],...
                    'LineWidth',lw,'Color',clr)
            end
        end
    end
end
axis image
if plot3D
    view([0 0]);
end

%% Add molecule labels
c = 0;
for g = 1:max(G)
    ig = find(G==g);
    c = c + 1; if c > size(cs,1); c = 1; end
    
    TextProperties =  {'FontSize',8,'HorizontalAlignment','center',...
        'Color',cs(c,:),'BackgroundColor','w','EdgeColor','none',...
        'Margin',0.5,'LineWidth',0.5};
    
    for m = ig'
        if plot3D
            text([mg(m).x],[mg(m).y],[mg(m).z],names(m),TextProperties{:})
        else
            text([mg(m).x],[mg(m).z],names(m),TextProperties{:})
        end
    end
%     set(gca,'xdir','rev')

end

