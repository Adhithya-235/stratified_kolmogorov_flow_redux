function MOV = animate_slices(X,Z,PHI,clabel,time,nf)

%   Animates slices of the given data. 

%% INITIALIZE FIGURE

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])

%% PREALLOCATE FRAMES

MOV(nf) = struct('cdata',[],'colormap',[]);

%% COMPUTE COLOUR AXIS LIMITS

%Cmax = max(max(max(max(PHI))));
%Cmin = min(min(min(min(PHI))));
Cmax = 3;
Cmin = -3;

%% ANIMATE SLICES

for i=1:nf
   fprintf('Accessing frame %d of %d.\n',i,nf)
   help_plot_fields(X, Z, squeeze(PHI(:, :, i)), [Cmin,Cmax], clabel, time(i));
   drawnow
   MOV(i) = getframe(f); 
end

end