function d_qhat = diffn(qhat, K, varargin)
% Applies the derivative operator to the unwrapped 2D field qhat 
% and returns the unwrapped field d_qhat. The field qhat is given in 
% spectral space, and so is d_qhat. K is a wavenumber matrix. Enter Kx if
% an x-derivative is to be taken, and Kz for a z-derivative. If the field
% to be differentiated is a perturbation field with a nonzero Floquet
% exponent, vararagin must contain two inputs, L, the period of the base 
% state in the derivative direction, and a, the Floquet modifier under
% consideration. 

%% GET NO. OF GRID POINTS IN x AND z

Nx = size(K,2);
Nz = size(K,1);

%% WRAP INPUT VECTOR

Qhat = reshape(qhat, Nz, Nx); 

%% COMPUTE DERIVATIVE FOR BASIC STATE FIELD (EMPTY varargin)

if isempty(varargin)
    d_Qhat = 1i*K.*Qhat;
    d_qhat = d_Qhat(:);
end

%% COMPUTE DERIVATIVE FOR PERTURBATION FIELD (2-element varargin)

if ~isempty(varargin)
    
    % Check number of optional inputs. 

    if length(varargin) ~= 2
        error('Exactly 2 optional inputs (L and a) are required for perturbation fields.')
    end

    % Define optional inputs

    L = varargin{1}; % base state period
    a = varargin{2}; % floquet modifier
    k = 2*pi/L;      % fundamental wavenumber

    % Calculate shift matrix

    I = ones(Nz,Nx);

    % Calculate derivative

    d_Qhat = 1i*(K + k*a*I).*Qhat;
    d_qhat = d_Qhat(:);

end

end