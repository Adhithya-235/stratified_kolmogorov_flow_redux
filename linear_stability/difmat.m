function [z, D] = difmat(N, m)
%=========================================================================%
%   Construct fourier differentiation matrices for my specific problem
%   following Dr. Chini's code "fourdif.m"
%   
%   N -> # of modes (even),  m -> Order of differentiation (either 1 or 2)
%
%   z -> Fourier grid, D -> Diffn. matrix
%=========================================================================%

%% PRELIMINARIES

h = 2*pi/N;
z = ((0:N-1)*h)';

%% CONVENIENCE VARIABLES

ii = (1:N-1)';
n1 = floor((N-1)/2);
n2 = ceil((N-1)/2);

%% START, m = 1

if m == 1
    topc = cot(0.5*h*(1:n2)');
    col1 = [0; 0.5*((-1).^ii).*[topc; -flipud(topc(1:n1))]];
    row1 = -col1;

%% START, m = 2

elseif m == 2
    
    topc = (csc(0.5*h*(1:n2)')).^2;
    f1 = -((2*pi*pi + h*h)/(6*h*h));
    col1 = [f1; -0.5*((-1).^ii).*[topc; flipud(topc(1:n1))]];
    row1 = col1;

%% START, m > 2

else
   N1 = floor((N-1)/2);
   N2 = (-N/2)*rem(m+1,2)*ones(rem(N+1,2));
   mwave=1i*[(0:N1) N2 (-N1:-1)];
   col1=real(ifft((mwave.^m).*fft([1 zeros(1,N-1)])));
   if rem(m,2)==0
       row1=col1;                           
   else
       col1=[0 col1(2:N)]';
       row1=-col1;                         
   end
end

%% CONSTRUCT MATRIX

D = toeplitz(col1,row1);

end