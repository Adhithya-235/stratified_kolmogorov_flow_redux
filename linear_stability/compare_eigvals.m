clear
close all
clc

%% SOLUTION PARAMETERS

Fr     = 0.02;
Rb     = [0.5,1,10,50];
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

eigsolns = struct; eigsolns.grate={}; eigsolns.phspd={}; eigsolns.kx={};
eigsolns(Ne) = eigsolns;

for s = 1:Ne
    eigsolns(s) = load(fnames(s),'grate','phspd','kx');
end

%% PLOT DOMINANT ROOTS AT EACH kx

fs = 30;
lw = 3;

f1 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
subplot(122)
hold on 
for s = 1:Ne
    plot(eigsolns(s).kx,abs(eigsolns(s).phspd(1,:)),'-',...
        'LineWidth', lw, 'MarkerSize',3)
end
ylabel('$\mathcal{I}[\sigma]/k_x$','Interpreter','latex')
xlabel('$Fr\,k_x$','Interpreter','latex')
legend('$Re_b = 0.5$','$Re_b = 1$','$Re_b = 10$','$Re_b = 50$',...
    'interpreter','latex')
axis square
xlim([0,1.75])
ylim([0,2])
grid on
box on
set(gca,'FontSize',fs,'LineWidth',lw,'GridLineWidth',2,...
    'ticklabelinterpreter','latex')

subplot(121)
hold on 
for s = 1:Ne
    plot(eigsolns(s).kx,eigsolns(s).grate(1,:),'-',...
        'LineWidth', lw, 'MarkerSize',3)
end
ylabel('$Fr\, \mathcal{R}[\sigma]$','Interpreter','latex')
xlabel('$Fr\, k_x$','Interpreter','latex')
axis square
xlim([0,1.75])
ylim([-0.1,0.1])
grid on
box on
set(gca,'FontSize',fs,'LineWidth',lw,'GridLineWidth',2,...
    'ticklabelinterpreter','latex')

sgtitle(sprintf('Dominant Eigenvalues, Fr = %.3f, Pr = %.3f', Fr, Pr),...
    'interpreter', 'latex', 'fontsize', fs)
drawnow

%% SAVE PLOT

pfolder   = 'images\';
p2folder  = 'eignvalue_comparisons';
parstring = sprintf('Fr_%1.3f_Pr_%1.2f_Rb_varying_zoomed.png', Fr, Pr);
mkdir([pfolder,p2folder])
saveas(f1,parstring);