%=========================================================================%
%   Plot linear stability growth rate curve for the 2D-SSKF problem.
%=========================================================================%

clear
close all
clc

%% INVISCIDITY TOGGLE OFF

inviscid = 0;

%% ENTER PROBLEM PARAMETERS

H  = 4*pi/3;
N  = 128;
Fr = 0.02;
Rb = 1;
Pr = 1;
% k  = (0.34/3)*[1,6,12,19]; 
k  = linspace(0.001, 6, 250);
kx = k/Fr;
m  = 3;

%% GET z AND DIFF MATRICES

[~, D1] = difmat(N,1); 
[x, D2] = difmat(N,2); 
scale   = H/(2*pi);
z       = scale*x;
D1      = (1/(scale))*D1;
D2      = (1/(scale^2))*D2;

%% LOAD BASIC STATE IF DATA PRESENT

fprintf("Loading basic state data from file.\n")
load("linstab_basic.mat",'um_target','bm_target')
Ub  = um_target;
Bb  = bm_target;
Ub2 = D2*Ub;
Bb1 = D1*Bb;

%% SET UP SOLUTION ARRAYS

sigs = zeros(2*N, length(kx));
vort = zeros(N, 2*N, length(kx));
buoy = zeros(N, 2*N, length(kx));
strm = zeros(N, 2*N, length(kx));

%% SOLVE PROBLEM FOR EACH k

parfor i = 1:length(kx)
   fprintf('Now I am calculating k = %1.3f.\n',kx(i));
   [sigs(:,i), vort(:,:,i), buoy(:,:,i), strm(:,:,i), ~] = ...
       solve_2dsskf_gevp(Ub, Ub2, Bb1, N, H, Fr, Pr, Rb, kx(i), m,...
       inviscid); 
end

%% GET GROWTH RATE AND PHASE SPEED

grate = real(sigs);
phspd = imag(sigs)./kx;

%% TRANSFORMATIONS

kx    = Fr*kx;
grate = grate*Fr;

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
ylim([-4,4])
grid on
box on
set(gca,'FontSize',fs,'LineWidth',lw)

yyaxis right
hold on 
for i = 1:nr
    plot(kx,grate(i,:),'.','MarkerSize',12)
end
ylabel('$Fr\, \mathcal{R}[\sigma]$','Interpreter','latex')
xlabel('$Fr\, k_x$','Interpreter','latex')
axis square
ylim([-1, 1])
grid on
box on
set(gca,'FontSize',fs,'LineWidth',lw)
drawnow

%% SAVE PLOT AND DATA

parstring = sprintf('Fr_%1.3f_Reb_%3.2f_Pr_%1.2f', Fr, Rb, Pr);
filename = sprintf('solution_%s.mat',parstring);
fig      = sprintf('eigenvals.png');
fig2     = sprintf('eigenvals.fig');
saveas(f1,fig);
saveas(f1,fig2);
save(filename,"grate","phspd", "vort","buoy","strm","z",...
    "kx","Rb","Pr","Fr","m");


