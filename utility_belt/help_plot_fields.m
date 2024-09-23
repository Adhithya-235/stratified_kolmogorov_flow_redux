function [] = help_plot_fields(X,Z,PHI,clim,clabel,time)

%   Plot slices of PHI. Inputs are self-explanatory (heheh), 
%   apart from clim, which is to be entered as a two-element 
%   increasing vector. 

hold off
pcolor(X,Z,PHI)  
caxis manual
caxis(clim)
colormap(slanCM('seismic'))
shading interp
c1 = colorbar;
xlabel('$x$', 'interpreter', 'latex')
ylabel('$z$', 'interpreter', 'latex')
title(sprintf('Time = %.2f', time))
c1.Label.Interpreter = 'latex';
c1.Label.String = clabel;
c1.Limits = clim;
c1.FontSize = 20;
c1.Location = 'eastoutside';
hold on
set(gca, 'fontsize', 30)
axis tight
box on
set(gca, 'boxstyle', 'full', 'linewidth', 4)
daspect([2 1 1])

end