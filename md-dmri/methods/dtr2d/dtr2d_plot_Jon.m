function dtr2d_plot(S, xps, axh, axh2)
% function dtd_3d_plot(S, xps, axh, axh2)

if (nargin < 4), axh2 = []; end

ms = 5;
fs = 10;
lw = 1;

lw_contour = 1;
lw_axes = 2;
lw_plot = .5;
fs_axes = 15;

opt = mdm_opt();
opt = dtr2d_opt(opt);

opt.dtr2d.dmin = 1e-11;
opt.dtr2d.dmax = 3.5e-9;
opt.dtr2d.r2min = 1;
opt.dtr2d.r2max = 500;

% opt.dtd_3d.dmin = .1/max(xps.b);
% opt.dtd_3d.r2min = .1/max(xps.te);
% opt.dtd_3d.r2max = 1/min(xps.te);

m = feval('dtr2d_1d_data2fit', S, xps, opt);
S_fit = feval('dtr2d_1d_fit2data', m, xps)';


% dmin = opt.dtd_3d.dmin;
% dmax = opt.dtd_3d.dmax;
% ratiomin = dmin/dmax;
% ratiomax = dmax/dmin;
% r2min = opt.dtd_3d.r2min;
% r2max = opt.dtd_3d.r2max;

% 
% xmin = log10(dmin);
% xmax = log10(dmax);
% ymin = log10(1/ratiomax);
% ymax = log10(ratiomax);
% zmin = log10(r2min);
% zmax = log10(r2max);

xmin = -11;
xmax = -8;
ymin = -3;
ymax = 3;
zmin = -0.5;
zmax = 3;

[n,par,perp,theta,phi,r2,w] = dtr2d_dist2par(dtr2d_m2dtr2d(m));
s0 = sum(w);

cla(axh);
hold(axh, 'off');
hsc1 = plot(axh,1:xps.n,S,'o',1:xps.n,S_fit,'x');
set(hsc1,'MarkerSize',ms,'LineWidth',lw)
axis(axh,'tight')
set(axh,'XLim',xps.n*[-.1 1.1], 'YLim',s0*[-.1 1.1],...
    'Box','off','TickDir','out','TickLength',.02*[1 1],...
'FontSize',fs,'LineWidth',lw)
xlabel(axh,'Acq number','FontSize',fs)
ylabel(axh,'Signal','FontSize',fs)

cla(axh2);
hold(axh2, 'on');

if n>0
    xcos = cos(phi).*sin(theta);
    ycos = sin(phi).*sin(theta);
    zcos = cos(theta);

    iso = tm_eigvals2iso([par perp perp]);
    fa = tm_eigvals2fa([par perp perp]);
    
%     delta = dtr2d_eigenvals2delta(iso, [par perp perp]); %Calculates D_delta

    c.x = log10(iso);
%     c.x = delta;
    c.y = log10(par./perp);
    c.z = log10(r2);
    c.ms = 50*ms*sqrt(w/s0);
    c.bright = fa;
    c.r = abs(xcos);
    c.g = abs(ycos);
    c.b = abs(zcos);
    
%    col = [c.r c.g c.b];
    
    x = c.x(:); x(~isfinite(x)) = 0;
    y = c.y(:); y(~isfinite(y)) = 0;
    z = c.z(:); z(~isfinite(z)) = 0;
    a = w(:);
    
    for nc = 1:n
        col = c.bright(nc)*[c.r(nc) c.g(nc) c.b(nc)];
        h1 = scatter3(axh2,c.x(nc),c.y(nc),c.z(nc));
%         h1.MarkerFaceColor = col;
        h1.MarkerEdgeColor = col;
        h1.SizeData = c.ms(nc);
        hold on
        view(axh2, 45, 20);
        grid(axh2, 'on');
        
    end

    
    c.n = n; 
    c.x = x;
    c.y = y;
    c.w = a;

    dist_s.x = linspace(xmin, xmax, n*4);
    dist_s.y = linspace(ymin, ymax, n*4);
    dist_s.xsigma = 2*abs(dist_s.x(2) - dist_s.x(1));
    dist_s.ysigma = 2*abs(dist_s.y(2) - dist_s.y(1));

    dist_s = dist_2d_discrete2smooth(c, dist_s);

    Cont = contourc(dist_s.x, dist_s.y, dist_s.w', n);
    hcontour = [];
    count = 1;
    while count < length(Cont)
        numxy = Cont(2,count);
        xtemp = Cont(1,count+(1:numxy));
        ytemp = Cont(2,count+(1:numxy));
        h = plot3(axh2, xtemp,ytemp,zmin*ones(size(xtemp)),'k-','LineWidth',lw_contour);
        hcontour = [hcontour; h];
        count = count + numxy + 1;
    end
    
    c.n = n; %numel(c.x);
    c.x = x;
    c.y = z;
    c.w = a;

    dist_s.x = linspace(xmin, xmax, n*4);
    dist_s.y = linspace(zmin, zmax, n*4);
    dist_s.xsigma = 2*abs(dist_s.x(2) - dist_s.x(1));
    dist_s.ysigma = 2*abs(dist_s.y(2) - dist_s.y(1));

    dist_s = dist_2d_discrete2smooth(c, dist_s);

    Cont = contourc(dist_s.x, dist_s.y, dist_s.w', n);
    count = 1;

    while count < length(Cont)
        numxy = Cont(2,count);
        xtemp = Cont(1,count+(1:numxy));
        ytemp = Cont(2,count+(1:numxy));
        h = plot3(axh2, xtemp,ymax*ones(size(xtemp)),ytemp,'k-','LineWidth',lw_contour);
        hcontour = [hcontour; h];
        count = count + numxy + 1;
    end

    c.n = n; %numel(c.x);
    c.x = y;
    c.y = z;
    c.w = a;

    dist_s.x = linspace(ymin, ymax, n*4);
    dist_s.y = linspace(zmin, zmax, n*4);
    dist_s.xsigma = 2*abs(dist_s.x(2) - dist_s.x(1));
    dist_s.ysigma = 2*abs(dist_s.y(2) - dist_s.y(1));

    dist_s = dist_2d_discrete2smooth(c, dist_s);

    Cont = contourc(dist_s.x, dist_s.y, dist_s.w', n);
    count = 1;

    while count < length(Cont)
        numxy = Cont(2,count);
        xtemp = Cont(1,count+(1:numxy));
        ytemp = Cont(2,count+(1:numxy));
        h = plot3(axh2, xmin*ones(size(xtemp)),xtemp,ytemp,'k-','LineWidth',lw_contour);
        hcontour = [hcontour; h];
        count = count + numxy + 1;
    end
    
    
end

% 
% set(axh2,'XLim',[xmin xmax], 'YLim',[-3 3],'ZLim',[zmin zmax],'YAxisLocation','right',...
% 'XTick',[xmin:1:xmax],'YTick',[-3:1:3],'ZTick',[zmin:.25:zmax],'TickDir','out','TickLength',.02*[1 1],...
% 'FontSize',fs,'LineWidth',lw,'Box','on')
% 

set(axh2,'XLim',[xmin xmax], 'YLim',[ymin ymax],'ZLim',[zmin zmax],'YAxisLocation','right',...
'XTick',[xmin:1:xmax],'YTick',[ymin:1:ymax],'ZTick',[zmin:.5:zmax],'TickDir','out','TickLength',.02*[1 1],...
'FontSize',fs,'LineWidth',lw,'Box','on')
set(hcontour,'Color',.5*[1 1 1])
axis(axh2,'square')
xlabel(axh2,'size, log(\itD\rm_{iso}/m^2s^-^1)','FontSize',fs, 'Rotation',-20)
ylabel(axh2,'                     shape, log(\itD\rm_{||}/\itD\rm_{\perp})','FontSize',fs, 'Rotation',25)
zlabel(axh2,'log(R_{2}/s^{-1})','FontSize',fs)
title(axh2,['orientation, [RGB]=[xyz]'],'FontSize',fs,'FontWeight','normal')


%       hsc1 = @(ax,x,y) scatter(ax,x,y,c.ms,col,'Marker','*','LineWidth',.03);
%       hsc2 = @(ax,x,y) scatter(ax,x,y,c.ms,col,'filled','Marker','o','LineWidth',.01);
%       
%       [axh3,~,~] = plotyy(axh2,c.x,c.z,c.x,c.y,hsc1,hsc2);  % draw markers with handle
% %       set(hAx,'nextplot','add')                 % "hold on" for multiple axes
% 
% set([axh3(1); axh3(2)],'XLim',[xmin xmax],'XTick',[-11:.5:-8],'TickDir','out',...
%     'TickLength',.02*[1 1],'FontSize',fs,'LineWidth',lw)
% 
% set(axh2,'XLim',[xmin xmax],'XTick',-11:.5:-8,'TickDir','out',...
%     'TickLength',.02*[1 1],'FontSize',fs,'LineWidth',lw,'Box','on')
% 
% set(axh3(2),'YLim',[ymin ymax],'YTick',-2:.5:2,'YColor','b')
% set(axh3(1),'YLim',[-.5 1.5],'YTick',-.5:.5:1.5,'YColor','r')
% 
% 
% 
% axis([axh2; axh3(1); axh3(2)],'square')
% xlabel(axh2,'size, log(\itD\rm_{iso} / m^2s^-^1)','FontSize',fs)
% ylabel(axh3(1),'log(\it{R}\rm_2 / s^-^1)','FontSize',fs)
% ylabel(axh3(2),'shape, log(\itD\rm_{||} / \itD\rm_{\perp})','FontSize',fs)
% 
% 
% %         h2 = plot(axh2,c.x(nc),c.z(nc),'x','LineWidth',.01);
% % plot(axh2,[log10(dmin) log10(dmax)],log10(ratiomax)*[1 1],'k-','LineWidth',lw);
% % 
% % axis(axh2,'square')
% % xlabel(axh2,'size, log(\itD\rm_{iso} / m^2s^-^1)','FontSize',fs)
% % ylabel(axh2,'shape, log(\itD\rm_{||} / \itD\rm_{\perp})    log(\it{R}\rm_2 / s^-^1)','FontSize',fs)
% title(axh2,['orientation, [RGB]=[xyz]'],'FontSize',fs,'FontWeight','normal')
