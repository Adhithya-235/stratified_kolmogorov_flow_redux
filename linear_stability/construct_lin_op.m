function [L, B, z] = construct_lin_op(N, H, Fr, Pr, Rb, k, Ub, Ub2, Bb1, inviscid)

%=========================================================================%
% Construct discretized linear linear operators to solve the generalized
% eigenvalue problem LV = SBV, where S is a growth rate. L is given by the
% block matrix [He -Ik Jp; 0 Hb Kp; -I 0 Hp]. All blocks are N x N. B is
% given by [I 0 0; 0 I 0; 0 0 0]. The state variable is [vort,buoy,strm].
% The variable 'inviscid' toggles inviscidity. Set it to 1 to take Re, Pe
% to infinity.
%=========================================================================%

%% GRID AND DIFFERENTIATION MATRICES

[~, D1] = difmat(N, 1);
[x, D2] = difmat(N, 2);

%% SCALE GRID AND DMs

scale = H/(2*pi);

z  = scale*x;
D1 = (1/scale)*D1;
D2 = (1/(scale^2))*D2;

%% CONSTRUCT BLOCK MATRICES - PRELIMINARIES

U   = diag(Ub);
D2U = diag(Ub2);
D1B = diag(Bb1);
I   = eye(N);
F2  = ((Fr*k)^2)*I;
Zm  = zeros(N);

%% CONSTRUCT BLOCKS in L 

He = (1/Rb)*(D2 - F2) - 1i*k*U;
Ik = 1i*k*eye(N);
Jp = 1i*k*D2U;
Hb = (1/(Rb*Pr))*(D2 - F2) - 1i*k*U;
Kp = 1i*k*(I + D1B);
Hp = D2 - F2;

%% OVERWRITE IF WE WANT INVISCID ROOTS

if inviscid==1
    He = - 1i*k*U;
    Hb = - 1i*k*U;
end

%% ASSEMBLE L and B

L = [He -Ik Jp; Zm Hb Kp; -I Zm Hp];
B = [I Zm Zm; Zm I Zm; Zm Zm Zm];

end