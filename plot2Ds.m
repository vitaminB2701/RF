function h = plot2Ds(data, Tw, varargin)
% PLOT2DDAS Contour plot of 2D spectra
%
% Synthax:
%   plot2Ds(Data,Tw)
%   plot2Ds(Data,Tw,Name,Value)
%   h = plot2Ds(___)
%
% Description
%   plot2Ds(Data,Tw) plots the 2D spectra in Data at delay times Tw.
%
%   Data is a struct with fields:
%       X   - X frequency axis (t3)
%       Y   - Y frequency axis (t1)
%       T   - delay time axis (Tw)
%       Abs - 3D array (Y,X,T) signal
%
%   h = plot2Ds(___) returns an array of figure handles.
%
% Name-Value Pair Arguments
%   PlotType      - plot type   ('contour', 'contourf', 'image',
%                                           'overlay', 'smooth')
%   Colormap      - colormap                ('parula, 'jet', ...)
%   FixedColor    - use a fixed color scale for all spectra (t|f)
%   FixedContour  - use fixed contour levels                (t|f)
%   CLevels       - number of contour levels             (scalar)
%   XLim          - X axis limits                       [min max]
%   YLim          - Y axis limits                       [min max]
%   Unit          - frequency axis units            ('nm'|'cm-1')
%   Style         - figure style (see figstyle)          (string)
%   FigureSize    - figure width and height in cm  [width height]
%   FontSize      - axes font size                       (scalar)
%   LineWidth     - axes line width                      (scalar)
%   XLabel        - X label                              (string)
%   YLabel        - Y label                              (string)
%   Title         - figure title                         (string)
%   Colorbar      - show colorbar                    ('on','off')
%   DiagColor     - diagonal color                    (or 'none')

%% Parse input parameters
p = inputParser;

p.addOptional('n',[]);
p.addParameter('FixedColor',true);
p.addParameter('FixedContour',false);
p.addParameter('CLevels',8);
p.addParameter('XLim',[]);
p.addParameter('YLim',[]);
p.addParameter('Unit','nm');
p.addParameter('Format','');
p.addParameter('Style','');
p.addParameter('FigureSize',[]);
p.addParameter('XLabel', 0);
p.addParameter('YLabel', 0);
p.addParameter('Title', 0);
p.addParameter('Colorbar','on');
p.addParameter('FontSize',[]);
p.addParameter('LineWidth',[]);
p.addParameter('PlotType','smooth');
p.addParameter('Colormap', '');
p.addParameter('DiagColor', 'k');

p.parse(varargin{:});
pars = p.Results;

%% Validate input
if ~isstruct(data) || ~all(isfield(data,{'X','Y','T','Abs'}))
    error('Invalid input - Data must be a struct with fields X, Y, T, Abs.')
end
if ~isnumeric(Tw)
    error('Invalid input - Tw must be a number or numeric array.')
end

%% Initialize

% Set X,Y labels
if strcmp(pars.Unit,'nm')
    YUnit = '\lambda_\tau (nm)';
    XUnit = '\lambda_t (nm)';
elseif  strcmp(pars.Unit,'cm')
    YUnit = '\omega_\tau (cm^-^1)';
    XUnit = '\omega_t (cm^-^1)';
elseif  strcmp(pars.Unit,'THz')
    YUnit = '\omega_\tau (THz)';
    XUnit = '\omega_t (THz)';
elseif  strcmp(pars.Unit,'eV')
    YUnit = '\omega_\tau (eV)';
    XUnit = '\omega_t (eV)';    
else
    YUnit = strcat('\omega_\tau (', pars.Unit, ')');
    XUnit = strcat('\omega_t (', pars.Unit, ')');
end

if ischar(pars.XLabel)
    XUnit = XLabel;
end
if ischar(pars.YLabel)
    YUnit = YLabel;
end

% Select delay times
Ti = zeros(1,length(Tw));
for i = 1:length(Tw)
    [~,Ti(i)] = min(abs(data.T - Tw(i)));
end
Ti = unique(Ti);
if length(Ti) < 1
    error 'Could not find the delay times Tw in the data set.'
end
Abs = real(data.Abs(:,:,Ti));

% Calculate fixed color scheme
if pars.FixedContour
    pars.FixedColor = true;
end

if pars.FixedColor
    % find amplitude range
    dasmin = min(Abs(:));
    dasmax = max(Abs(:));
    clst = 0; % contour level step
    if pars.FixedContour
        clst = (dasmax - dasmin) / (pars.CLevels+1);
        pars.CLevels = dasmin:clst:dasmax-clst;
    end
end

%% Plot
handles = [];
for i = 1:numel(Ti)
    h = figure();
    handles = [handles, h];
    
    % Set figure format
    if ~isempty(pars.Format)
        formatfig(pars.Format)
    end
    if ~isempty(pars.Style)
        setfigstyle(h,pars.Style)
    end
    
    ax = gca;
    if isempty(pars.LineWidth)
        LineWidth = get(gcf,'DefaultAxesLineWidth');
    else
        LineWidth = pars.LineWidth;
        set(gcf,'DefaultAxesLineWidth',LineWidth);
        set(gcf,'DefaultLineLineWidth',LineWidth);
    end
    if isempty(pars.FontSize)
        FontSize = get(gcf,'DefaultAxesFontSize');
    else
        FontSize = pars.FontSize;
        set(gcf,'DefaultAxesFontSize',FontSize);
    end
    
    if ~isempty(pars.FigureSize)
        setfigsize(pars.FigureSize,'centimeters');
    end
    
    % Create contour plot
    X = data.X; Y = data.Y; Z = Abs(:,:,i);
    switch lower(pars.PlotType)
        case 'contourf'
            contourf(ax,X,Y,Z,pars.CLevels,'LineWidth',LineWidth);
        case 'contour'
            contour(ax,X,Y,Z,pars.CLevels,'LineWidth',LineWidth);
        case 'image'
            [Xi,Yi,Zi] = interpolate(X,Y,Z);
            imagesc(ax,Xi,Yi,Zi);
            shading interp
        case 'pcolor'
            pcolor(ax,X,Y,Z)
            shading interp
            set(gca,'layer','top')
        case 'overlay'
            [Xi,Yi,Zi] = interpolate(X,Y,Z);
            imagesc(ax,Xi,Yi,Zi);
            hold on
            contour(ax,X,Y,Z,pars.CLevels,'LineWidth',LineWidth, 'Color', 'k');
        case 'smooth'
            [Xi,Yi,Zi] = interpolate(X,Y,Z);
            imagesc(ax,Xi,Yi,Zi);
            hold on
            
            % Smooth data for contour plot
            smZ = Z;
            for k = 1:numel(Y)
                smZ(k,:) = smooth(Z(k,:),5);
            end
            contour(ax,X,Y,smZ,pars.CLevels,'LineWidth',LineWidth, 'Color', 'k');
            
        otherwise
            error('Invalid argument. Type "help plot2Ds" for a list of valid plot types.')
    end
    
    axis image
    
    
    % Axes scale and gridlines
    if ~isempty(pars.XLim)
        xlim(pars.XLim)
    end
    if ~isempty(pars.YLim)
        ylim(pars.YLim)
    end
    xl = xlim; yl = ylim;
    dmin = max([xl(1) yl(1)]);
    dmax = min([xl(2) yl(2)]);
    
    % Diagonal line
    if ~strcmpi(pars.DiagColor,'none')
        line([dmin dmax], [dmin dmax], 'Color',pars.DiagColor, 'LineStyle','--','LineWidth',LineWidth);
    end
    
    grid off
    
    if strcmp(pars.Unit,'nm')
        ax.XDir = 'rev';
        ax.YDir = 'rev';
    else
        ax.YDir = 'normal';
    end
%     if strcmp(pars.Unit,'THz')
%         ax.YDir = 'normal';
%     end
    % Set color map
    if pars.FixedColor
        caxis([dasmin dasmax-clst]);
    end
%             caxis([-50 7]);
    
    if isempty(pars.Colormap)
        load colormap256
    elseif ischar(pars.Colormap)
        cmap = colormap(lower(pars.Colormap));
        cmap = cmap(end:-1:1,:);
    else
        cmap = colormap(pars.Colormap);
    end
    clim = caxis;
    numc = fix(0.5*size(cmap,1));
    if clim(1) >= 0
        cmap = cmap(numc:end,:);
    elseif clim(2) <= 0
        cmap = cmap(1:numc,:);
    elseif clim(2) > - clim(1)
        ii = 1+fix((1+clim(1)/clim(2))*numc);
        cmap = cmap(ii:end,:);
    else
        ii = round((1-clim(2)/clim(1))*numc);
        cmap = cmap(1:ii,:).^1;
    end
    colormap(cmap);
    
    % Axes title and labels
    xlabel(XUnit);
    ylabel(YUnit);
    
    if ischar(pars.Title)
        title(pars.Title)
    else
        if Tw(i) < 100
            title(sprintf('T_w = %1.4g ps', Tw(i)));
        elseif Tw(i) < 1000
            title(sprintf('T_w = %0.4g ps', Tw(i)));
        else
            title(sprintf('T_w = %0.4g ns',  Tw(i)/1000));
        end
    end
    
    % Color bar
    colorbar('LineWidth',LineWidth,'FontSize',FontSize,'TickLength',ax.TickLength(1))
    if strcmpi(pars.Colorbar,'off')
        axPos = ax.Position;
        colorbar hide
        ax.Position = axPos;
    end
end

if nargout > 0
    h = handles;
end
end

function [Xi,Yi,Zi] = interpolate(X,Y,Z)
% INTERPOLATE - 2D interpolation
n = 600;
[Xg,Yg] = meshgrid(X,Y);
xmin = min(X); xmax = max(X);
Xi = xmin:(xmax-xmin)/(n-1):xmax;
ymin = min(Y); ymax = max(Y);
Yi = ymin:(ymax-ymin)/(n-1):ymax;
[Xq,Yq] = meshgrid(Xi,Yi);
Zi = interp2(Xg, Yg, Z, Xq, Yq);
end