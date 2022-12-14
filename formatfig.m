% Format figure for presentation
%
% Syntax: formatfig(type, colors)
%
% type can be one of the following:
%   journal, colorjournal, poster, slide
%
% colors specifies the default color order as in colorscheme

% default type (if not specified) is colorjournal


function formatfig(style,colors)

% set default style
if nargin < 1
    style = 'default';
end

rootunits = get(0,'Units');
if ~strcmp(rootunits, 'pixels'), set(0,'Units','pixels'); end
dpi = get(0,'ScreenPixelsPerInch');
screensize = get(0,'ScreenSize');

f1 = gcf; a1 = gca;
switch lower(style)
    case 'default'
        set(f1,'Color','w');
        set(f1,'DefaultAxesLineWidth', 1);
        set(f1,'DefaultAxesTickLength', [0.02 0.02]);
        % set(f1,'DefaultAxesFontSize',16);
        box(gca,'on');    
        
    case 'journal'        
        set(f1,'DefaultAxesLineWidth', 1);
        set(f1,'DefaultLineLineWidth', 1);
        set(f1,'DefaultAxesTickLength', [0.02 0.02]);
        % set(f1,'DefaultAxesFontName','Arial');   
        set(f1,'DefaultAxesFontSize',16);   
        set(f1,'DefaultTextFontSize',16);
        set(f1,'DefaultLineMarkerSize',6);        
        set(f1,'DefaultAxesXColor','k','DefaultAxesYColor','k','DefaultAxesZColor','k');
%       set(f1,'DefaultAxesLabelFontSizeMultiplier',1);
        setfiguresize([13 10], 'centimeters')        
        box(gca,'on');
        set(gcf,'Color','w');
        set(a1,'color','none');
        
    case 'springer'        
        set(f1,'DefaultAxesLineWidth', 0.8);
        set(f1,'DefaultLineLineWidth', 1);
        set(f1,'DefaultAxesTickLength', [0.02 0.02]);
        set(f1,'DefaultAxesFontSize',10);   
        set(f1,'DefaultTextFontSize',10);
        set(f1,'DefaultLineMarkerSize',4);        
        set(f1,'DefaultAxesXColor','k','DefaultAxesYColor','k','DefaultAxesZColor','k');
        setfiguresize([9.1 7], 'centimeters')        
        set(f1,'Color','w');
        box(a1,'on');
        set(a1,'color','none');        
        
    case 'colorjournal'
        set(f1,'Color','w');
        set(f1,'DefaultAxesLineWidth', 1);
        set(f1,'DefaultLineLineWidth', 1.5);
        set(f1,'DefaultAxesTickLength', [0.02 0.02]);
        set(f1,'DefaultAxesFontSize',16);
        set(f1,'PaperUnits', 'centimeters', 'PaperType','A4');
        set(f1,'DefaultAxesXColor','k','DefaultAxesYColor','k','DefaultAxesZColor','k');
        set(f1,'PaperPosition',[2.5 8 16 12.8]);
        set(f1,'Position',[340 212 600 480]);              
        box(gca,'on');    
        
    case 'poster'
        set(f1,'Color','w');
        set(f1,'DefaultAxesLineWidth', 2);
        set(f1,'DefaultLineLineWidth', 2);
        set(f1,'DefaultAxesTickLength', [0.02 0.02]);
        set(f1,'DefaultAxesFontName','Arial');        
        set(f1,'DefaultAxesFontWeight','normal');
        set(f1,'DefaultAxesFontSize',24);
        setfiguresize([20 16],'centimeters')   
        box(gca,'on');

    case 'largeposter'
        set(f1,'Color','w');
        set(f1,'DefaultAxesLineWidth', 2);
        set(f1,'DefaultLineLineWidth', 2);
        set(f1,'DefaultAxesTickLength', [0.02 0.02]);
        set(f1,'DefaultAxesFontName','Arial');        
        set(f1,'DefaultAxesFontWeight','normal');
        set(f1,'DefaultAxesFontSize',16.5);
        set(f1,'DefaultAxesXColor','k','DefaultAxesYColor','k','DefaultAxesZColor','k');
        setfiguresize([12 9],'centimeters')
        box(gca,'on');

        case 'largeposter2'
        set(f1,'Color','w');
        set(f1,'DefaultAxesLineWidth', 2);
        set(f1,'DefaultLineLineWidth', 2);
        set(f1,'DefaultAxesTickLength', [0.02 0.02]);
        set(f1,'DefaultAxesFontName','Arial');        
        set(f1,'DefaultAxesFontWeight','normal');
        set(f1,'DefaultAxesFontSize',14.5);
        setfiguresize([14 11.2233],'centimeters')
        box(gca,'on');
        
     case 'journalsmall'
        set(f1,'DefaultAxesLineWidth', 0.7);
        set(f1,'DefaultLineLineWidth', 0.7);
        set(f1,'DefaultAxesTickLength', [0.02 0.02]);
        set(f1,'DefaultAxesFontSize',9);        
        set(f1,'DefaultLineMarkerSize',2.5);        
        set(f1,'DefaultAxesXColor','k','DefaultAxesYColor','k','DefaultAxesZColor','k');
        setfiguresize([9 7.2],'centimeters')
        box(gca,'on');
        set(a1,'color','none');
        
 case 'slide'
        set(f1,'DefaultAxesLineWidth', 1.5);
        set(f1,'DefaultLineLineWidth', 1.5);
        set(f1,'DefaultAxesTickLength', [0.02 0.02]);
        set(f1,'DefaultAxesFontSize',16);
        setfiguresize([640 480],'pixels')        
        box(gca,'on');
        set(a1,'color','none');

 case 'slide2d'
        set(f1,'DefaultAxesLineWidth', 0.5);
        set(f1,'DefaultLineLineWidth', 0.5);
        set(f1,'DefaultAxesTickLength', [0.02 0.02]);
        setfiguresize([11 8.5],'centimeters')        
        box(gca,'on');
        
    case 'dark'
        % Dark Background
        set(gcf,'Color','none');
        set(gca,'Color','none', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w', 'GridColor', 'w', 'MinorGridColor', 'w');       
       
    otherwise
        if nargin == 1
            set(f1,'DefaultAxesColorOrder',colorscheme(style));
        else
            error 'The supplied argument is not a valid formatting style'
        end
end
        if exist('colors','var')
            set(f1,'DefaultAxesColorOrder',colorscheme(colors));
        end
end

function setfiguresize(figuresize, units)
% sets the same figure size on screen and paper
% setfiguresize(figuresize,units)
% figuresize - desired figure size [width height]
% units - 
%   {inches} | centimeters | points

if ~exist('units','var')
    units = 'inches';
end

rootunits = get(0,'Units');
if ~strcmp(rootunits, 'pixels'), set(0,'Units','pixels'); end
dpi = get(0,'ScreenPixelsPerInch');
screensize = get(0,'ScreenSize');

% px - figure size in pixels
switch units
    case 'inches'
        px = figuresize * dpi;
    case 'centimeters'
        px = figuresize / 2.54 * dpi;
    case 'points'
        px = figuresize / 72 * dpi;
    case 'pixels'
        px = figuresize;
        figuresize = figuresize / dpi;
        units = 'inches';
    otherwise
        error('units can only be pixels, inches, centimeters or points');
end

set(gcf, 'Units', 'pixels', 'Position', [(screensize(3)-px(1))/2, (screensize(4)-px(2))/2, px], ...
         'PaperUnits', units, 'PaperSize', figuresize, 'PaperPosition', [0 0 figuresize]);
     
end