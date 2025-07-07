clear
close all
clc

%% PLOT AESTHETICS

fs   = 24;
lw   = 3;

%% PARAMETERS

Reb   = [1:20, 30:10:100];
Fr    = 0.01;
alpha = 0;
beta  = 0;

%% GET FILE NAMES

SolnFilePath = cell(length(Reb),1);
for fileidx = 1:length(Reb)
    ParentDir             = sprintf('../solutions_branch_1_Pr_1');
    subDir                = sprintf('Reb%.2f_alpha%.2f',Reb(fileidx), alpha);
    SolnFile              = sprintf('spectrum_Reb%.2f_alpha%.2f.mat',Reb(fileidx), alpha);
    SolnFilePath{fileidx} = sprintf('%s/%s/%s',ParentDir,subDir,SolnFile);
end

%% PREALLOCATE SOLUTION

nReb      = length(Reb);
nEigs     = 20;
list_eigs = NaN(nEigs, 1);

%% PREALLOCATE FRAMES

MOV(nReb) = struct('cdata',[],'colormap',[]);

%% PLOTTING LOOP

for idx = 1:nReb

    %% LOAD FILE AND EXTRACT EIGENVALUES

    data      = load(SolnFilePath{idx});
    list_eigs = data.eigvals(1:nEigs);

    %% PLOT EIGENVALUES

    espec = figure('Color', 'w', 'WindowState','maximized');
    hold off;
    scatter(real(list_eigs), imag(list_eigs), 40, 'filled');
    xlabel('$\mathrm{Fr} \cdot \mathrm{Re}(\sigma)$', 'Interpreter', 'latex', 'FontSize', 30);
    ylabel('$\mathrm{Fr} \cdot \mathrm{Im}(\sigma)$', 'Interpreter', 'latex', 'FontSize', 30);
    title(['$Re_b = $', num2str(Reb(idx))], 'interpreter', 'latex');
    xlim([-0.3, 0.3])
    ylim([-1, 1])
    axis square;
    grid on;
    box on
    set(gca, 'FontSize', 30, 'LineWidth', 3, 'TickLabelInterpreter', 'latex')
    drawnow
    
    %% CAPTURE FRAME

    MOV(idx) = getframe(espec); 

end

%% CREATE AVI FILE

fprintf('Writing movie to file...\n')
filename = sprintf('eigenspectrum_alpha%.2f.avi',alpha);
v = VideoWriter(filename); 
v.FrameRate = 1;
open(v)
for i = 1:nReb
   fprintf('Writing frame %d of %d.\n',i,nReb)
   writeVideo(v, MOV(i)) 
end
close(v)



