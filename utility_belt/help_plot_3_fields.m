function help_plot_3_fields(x, z, field1, field2, field3, label1, label2, label3, cmap, fs, lw)

%% PLOT FIELD 1 

subplot(3,1,1);
pcolor(x, z, field1); 
colormap(cmap)
shading interp;
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = label1;
cb.Label.FontSize = fs;
cb.FontSize = fs - 4; 
cb.LineWidth = lw / 2; 
cb.TickLabelInterpreter = 'latex';
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', fs-2);
ax1 = gca;
ax1.FontSize = fs - 4;
ax1.LineWidth = lw/2; 
ax1.TickLabelInterpreter = 'latex';
daspect([1 1 1])

%% PLOT FIELD 2

subplot(3,1,2);
pcolor(x, z, field2); 
colormap(cmap)
shading interp;
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = label2;
cb.Label.FontSize = fs;
cb.FontSize = fs - 4; 
cb.LineWidth = lw / 2; 
cb.TickLabelInterpreter = 'latex';
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', fs-2);
ax1 = gca;
ax1.FontSize = fs - 4;
ax1.LineWidth = lw/2; 
ax1.TickLabelInterpreter = 'latex';
daspect([1 1 1])

%% PLOT FIELD 3

subplot(3,1,3);
pcolor(x, z, field3); 
colormap(cmap)
shading interp;
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = label3;
cb.Label.FontSize = fs;
cb.FontSize = fs - 4; 
cb.LineWidth = lw / 2; 
cb.TickLabelInterpreter = 'latex';
xlabel('$x/Fr$', 'Interpreter', 'latex', 'FontSize', fs-2);
ylabel('$z$', 'Interpreter', 'latex', 'FontSize', fs-2);
ax1 = gca;
ax1.FontSize = fs - 4;
ax1.LineWidth = lw/2; 
ax1.TickLabelInterpreter = 'latex';
daspect([1 1 1])

end