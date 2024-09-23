function [datam, datap] = get_pert_fields(data,x,Lx,unwrap)

%   Performs Reynolds decomposition of data by averaging in the streamwise direction. 
%   The trapezoidal rule is used, which is exact for periodic data provided the endpoints 
%   are included. If unwrap is set to 1, the fluctuation field loses its last column. This
%   is useful for Fourier analysis of the perturbation field. 

%% CALCULATE STREAMWISE MEAN

datam = trapz(x, data, 2)/Lx;                    % Streamwise average

%% CALCULATE FLUCTUATION

datap = data - datam;

%% UNWRAP

if unwrap == 1
   datap(:, end, :) = [];
   datap(end, :, :) = [];
end

end