

%% DNS SIMULATION PARAMETERS

if alpha == 0
    alpha = 1;
end
Lx = alpha*0.03;
Lz = 2*pi/3;

%% READ DNS DATA

[x, z, ~, ~]              = get_space_data(folder_name, data_folder, file_name, 1);
[t, ~, ~, b, vort, nf]    = get_field_data(folder_name, data_folder, file_name, stride, svec, 0);

%% ECS FILE PARAMETERS

ecsFolder = sprintf('../secondary_stability/import');
ecsFile   = sprintf('ECS_real_field_for_Reb=%d.mat',Rb);
ecsPath   = sprintf('%s/%s',ecsFolder,ecsFile);

%% READ ECS DATA

ecs          = load(ecsPath);
ecs.omega    = interpft(interpft(repmat(ecs.omega, 1, 1/alpha), Nz, 1), Nx, 2);
ecs.Buoyancy = interpft(interpft(repmat(ecs.Buoyancy, 1, 1/alpha), Nz, 1), Nx, 2); 

%% CALCULATE PERTURBATIONS AND PERTURBATION ENERGIES

vortp = vort - ecs.omega;
bp    = b - ecs.Buoyancy;
enstr = vortp.^2;
poten = bp.^2;

%% WRAP DATA BEFORE VOLUME AVERAGE

poten = cat(2, poten, poten(:, 1, :));
enstr = cat(2, enstr, enstr(:, 1, :));
poten = cat(1, poten, poten(1, :, :));
enstr = cat(1, enstr, enstr(1, :, :));
Lx    = x(end);
Lz    = z(end);

%% VOLUME AVERAGE PERTURBATION ENERGIES

disp('Starting volume average.')
genstr = calc_volm_avg(enstr,x,Lx,z,Lz);
disp('Done with enstrophy.')
gpoten = calc_volm_avg(poten,x,Lx,z,Lz);
disp('Ending volume average.')

%% PLOT TIMESERIES

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
hold on
plot(t, genstr, '-o', 'linewidth', 3)
plot(t, gpoten, '-o', 'linewidth', 3)
xlabel('$t$', 'interpreter', 'latex')
ylabel('Pert Energy', 'interpreter', 'latex')
legend('Enstrophy','PE', 'SSGR', 'interpreter', 'latex')
set(gca, 'fontsize', 30)
xlim([t(1), t(end)])
grid on
box on
set(gca, 'linewidth', 5, 'XScale', 'linear', 'YScale', 'log')
drawnow
    
%% SAVE PLOT

saveas(f, sprintf('../%s/%s/plots/timeseries/pertenergy_timeseries.fig', folder_name, data_folder)) 
saveas(f, sprintf('../%s/%s/plots/timeseries/pertenergy_timeseries.png', folder_name, data_folder)) 
