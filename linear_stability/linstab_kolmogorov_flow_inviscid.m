%=========================================================================%
%   Plot linear stability growth rate curve for the 2D-SSKF problem.
%   Inviscid case. Ignores Rb and Pr but essential to specify them.
%=========================================================================%

clear
close all
clc

%% INVISCIDITY TOGGLE OFF

inviscid = 1;

%% ENTER PROBLEM PARAMETERS

H  = 4*pi/3;
N  = 256;
Fr = 0.02;
Rb = 1;
Pr = 1;
kx = linspace(1000, 2000, 20);
m  = 3;

%% SET UP SOLUTION ARRAYS

sigs = zeros(2*N, length(kx));
vort = zeros(N, 2*N, length(kx));
buoy = zeros(N, 2*N, length(kx));
strm = zeros(N, 2*N, length(kx));

%% SOLVE PROBLEM FOR EACH k

parfor i = 1:length(kx)
   fprintf('Now I am calculating k = %1.3f.\n',kx(i));
   [sigs(:,i), vort(:,:,i), buoy(:,:,i), strm(:,:,i), z] = ...
       solve_2dsskf_gevp(N, H, Fr, Pr, Rb, kx(i), m, inviscid); 
end

%% GET GROWTH RATE AND PHASE SPEED

grate = real(sigs);
phspd = imag(sigs)./kx;

%% PLOT EIGENVALUES IN COMPLEX PLANE

yl = 4;
fs = 20;
lw = 3;
nr = 1;

f1 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
yyaxis left
hold on 
for i = 1:nr
    plot(kx,phspd(i,:),'.','MarkerSize',12)
end
ylabel('$\mathcal{I}[\sigma]/k_x$','Interpreter','latex')
xlabel('$k_x$','Interpreter','latex')
axis square
ylim([-1,1])
grid on
box on
set(gca,'FontSize',fs,'LineWidth',lw,'GridLineWidth',2)

yyaxis right
hold on 
for i = 1:nr
    plot(kx,grate(i,:),'.','MarkerSize',12)
end
ylabel('$\mathcal{R}[\sigma]$','Interpreter','latex')
xlabel('$k_x$','Interpreter','latex')
axis square
ylim([-25,25])
grid on
box on
set(gca,'FontSize',fs,'LineWidth',lw,'GridLineWidth',2)