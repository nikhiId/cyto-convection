clc; clear all;

% nucleus radius (must be positive!)
kappa = 0.43;

% nucleus eccentricity
% for best results keep |e| > 0.01 (since bi-spherical coordinates become singular in the concentric limit)
e = -0.3;

% number of bi-spherical harmonic modes
% increase N_BiSp until desired convergence is achieved
N_BiSp = 15;
% too high N_BiSp can cause ill-conditioned matrices

% number of xi-points over which to discretize the functions U_n(xi)
% increase M_xi until desired convergence is achieved
M_xi = 100;

% velocity component to be plotted on the y = 0 plane; can be 'x' or 'z' or
% 'mag' for the magnitude of the velocity
plot_vel = 'mag';

% function returns:
% 1. the minimum (uz_min) and maximum (uz_max) vertical velocity inside the cell
% 2. the bi-spherical eigenfunctions/velocity modes U_n(xi)
% function plots the velocity component 'plot_vel' (entered in line-21) on the y=0 plane
[uz_max, uz_min, Un] = axisymm_flow_solve_BiSp(e, kappa, N_BiSp, M_xi, plot_vel);