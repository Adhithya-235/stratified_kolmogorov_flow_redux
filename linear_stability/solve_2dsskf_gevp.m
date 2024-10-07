function [sigm, vort, buoy, strm, z] = solve_2dsskf_gevp(Ub, Ub2, Bb1, N, H, Fr, Pr, Rb, k, m, inviscid)

%   Set up and solve a generalized eigenvalue problem for 2D Strongly
%   Stratified Kolmogorov Flow. 

%% CONSTRUCT LINEAR OPERATORS

[L, B, z] = construct_lin_op(N, H, Fr, Pr, Rb, k, Ub, Ub2, Bb1, inviscid);

%% SOLVE GEVP for MAX GROWTH RATE

[sigm, vort, buoy, strm] = find_maxgrowth(L, B);

end