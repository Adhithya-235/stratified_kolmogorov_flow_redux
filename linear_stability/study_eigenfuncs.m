clear
close all
clc

%% SOLUTION PARAMETERS

Fr     = 0.02;
Rb     = 50;
Pr     = 1;
Ne     = length(Rb);

%% SAVE FILENAMES

parent = 'solutions';
subfld = 'parameter_sweep';
direct = string.empty;
fnames = string.empty;
object = string.empty;

for s = 1:Ne
    direct(s) = sprintf('Fr_%1.3f_Reb_%3.2f_Pr_%1.2f', Fr, Rb(s), Pr);
    object(s) = sprintf('solution_Fr_%1.3f_Reb_%3.2f_Pr_%1.2f.mat', Fr, Rb(s), Pr);
    fnames(s) = sprintf('%s/%s/%s/%s',parent,subfld,direct(s),object(s));
end

%% LOAD FILES

eigsolns = struct; 
eigsolns.grate={}; 
eigsolns.phspd={}; 
eigsolns.kx={};
eigsolns.buoy={};
eigsolns.z={};
eigsolns(Ne) = eigsolns;

for s = 1:Ne
    eigsolns(s) = load(fnames(s),'grate','phspd','kx','buoy','z');
end

%% SET UP PLOT

fs = 30;
lw = 3;
indices = [35,58,126,177];
px      = eigsolns.kx(indices);
py      = eigsolns.grate(1,indices);
f1 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
tiledlayout('flow')

%%  GROWTH RATE PLOT

nexttile([1,4])
hold on 
for s = 1:Ne
    plot(eigsolns(s).kx,eigsolns(s).grate(1,:),'-',...
        'LineWidth', lw, 'MarkerSize',3)
end
plot(px,py,'k.','MarkerSize',25)
ylabel('$Fr\, \mathcal{R}[\sigma]$','Interpreter','latex')
xlabel('$Fr\, k_x$','Interpreter','latex')
ylim([-1,1])
grid on
box on
set(gca,'FontSize',fs,'LineWidth',lw,'GridLineWidth',2,...
    'ticklabelinterpreter','latex')
title(sprintf('Eigenfunctions, Fr = %.3f, Pr = %.3f, Rb = %.3f', Fr, Pr, Rb),...
    'interpreter', 'latex', 'fontsize', fs)
drawnow

%% PLOT BRANCH I EIGENFUNCTION

nexttile([2,1])
hold on 
for s = 1:Ne
    plot(abs(eigsolns(s).buoy(:,1,35)),eigsolns(s).z,'-',...
        'LineWidth', lw, 'MarkerSize',3)
end
ylabel('$z$','Interpreter','latex')
xlabel('$\left|\hat{b}_{I}\right|$','Interpreter','latex')
grid on
box on
axis tight
set(gca,'FontSize',fs,'LineWidth',lw,'GridLineWidth',2,...
    'ticklabelinterpreter','latex')
drawnow

%% PLOT BRANCH II EIGENFUNCTION

nexttile([2,1])
hold on 
for s = 1:Ne
    plot(abs(eigsolns(s).buoy(:,1,58)),eigsolns(s).z,'-',...
        'LineWidth', lw, 'MarkerSize',3)
end
xlabel('$\left|\hat{b}_{II}\right|$','Interpreter','latex')
grid on
box on
axis tight
set(gca,'FontSize',fs,'LineWidth',lw,'GridLineWidth',2,...
    'ticklabelinterpreter','latex')
drawnow

%% PLOT BRANCH III EIGENFUNCTION

nexttile([2,1])
hold on 
for s = 1:Ne
    plot(abs(eigsolns(s).buoy(:,1,126)),eigsolns(s).z,'-',...
        'LineWidth', lw, 'MarkerSize',3)
end
xlabel('$\left|\hat{b}_{III}\right|$','Interpreter','latex')
grid on
box on
axis tight
set(gca,'FontSize',fs,'LineWidth',lw,'GridLineWidth',2,...
    'ticklabelinterpreter','latex')
drawnow

%% PLOT BRANCH IV EIGENFUNCTION

nexttile([2,1])
hold on 
for s = 1:Ne
    plot(abs(eigsolns(s).buoy(:,1,177)),eigsolns(s).z,'-',...
        'LineWidth', lw, 'MarkerSize',3)
end
xlabel('$\left|\hat{b}_{IV}\right|$','Interpreter','latex')
grid on
box on
axis tight
set(gca,'FontSize',fs,'LineWidth',lw,'GridLineWidth',2,...
    'ticklabelinterpreter','latex')
drawnow

%% SAVE PLOT

pfolder   = 'images\';
p2folder  = 'eigenfunctions';
parstring = sprintf('Fr_%1.3f_Pr_%1.2f_Rb_%1.2f_buoy.png', Fr, Pr, Rb);
mkdir([pfolder,p2folder])
saveas(f1,parstring);