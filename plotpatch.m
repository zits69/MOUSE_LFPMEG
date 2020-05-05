function [linehandle patchhandle dothandle cfg] = plotpatch(y,timevec,patchcolor,cfg)
%[linehandle patchhandle dothandle cfg] = plotpatch(y,timevec,patchcolor,cfg)

nsubs = size(y,1);
ntime = size(y,2);
ncond = size(y,3);

if ndims(y) > 3
    error('plotpatch: Input Data Matrix has too many dimensions!');
end

if nargin < 4
   cfg           = [];
   cfg.linewidth = 4;
   cfg.alpha     = 0.5;
   cfg.plotdots  = false;
   cfg.plotline  = true;
end

if ~isfield(cfg,'linewidth')
    cfg.linewidth = 4;
end
if ~isfield(cfg,'alpha')
    cfg.alpha = 0.5;
end
if ~isfield(cfg,'plotdots')
    cfg.plotdots = false;
end
if ~isfield(cfg,'plotline')
    cfg.plotline = true;
end

if nargin < 3
    patchcolor = colormap(gray(ncond));
end

if size(patchcolor,1) == 1 & ncond > 1
   warning('plotpatch: Not enough colors for number of conditions! Plotting in grayscale.');
   patchcolor = colormap(gray(ncond));
end 

if nargin < 2
   timevec = 1:ntime; 
end

xp  = [timevec fliplr(timevec)];

patchhandle = [];
linehandle  = [];
dothandle   = [];
hold on
for icondition = 1:ncond
    ymu = squeeze(nanmean(y(:,:,icondition),1));
    ysd = squeeze(nanstd(y(:,:,icondition),[],1))./sqrt(size(y,1));    
    yp  = [ymu+ysd fliplr(ymu-ysd)];
    cp  = patchcolor(icondition,:);
    hold on
    patchhandle(icondition) = patch(xp,yp,cp,'EdgeColor','none','FaceAlpha',cfg.alpha);
    if cfg.plotline
        linehandle(icondition)  = plot(timevec,ymu,'k-','linewidth',cfg.linewidth,'color',cp);
    end
    if cfg.plotdots
       dothandle(icondition) = plot(timevec,ymu,'ko','Markerfacecolor',cp,'markeredgecolor',[0 0 0],'markersize',12,'linewidth',1);
    end
end

end