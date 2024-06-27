
% function takes as input the following:
% 1. nucleus radius (kappa) & eccentricity (e)
% 2. number of bi-spherical harmonic modes (N_modes)
% 3. number of xi-points over which to discretize the functions U_n(xi) (M)
% 4. velocity component to be plotted on the y = 0 plane (plot_vel)

% function returns the following:
% 1. the minimum (uz_min) and maximum (uz_max) vertical velocity inside the cell
% 2. the bi-spherical eigenfunctions U_n(xi)

% function plots the velocity component 'plot_vel' (entered in line-16) on the y=0 plane

function [uz_max, uz_min, Un] = axisymm_flow_solve_BiSp(e, kappa, N_modes, M, plot_vel)

if kappa + abs(e) < 1

    Ng = 200; N_chi_plot = 100;

    xi_i = -e/abs(e)*acosh( (1 - kappa^2 - abs(e)^2)/2/kappa/abs(e) );
    xi_o = -e/abs(e)*acosh( (1 - kappa^2 + abs(e)^2)/2/abs(e) );

    a = abs( sinh(xi_o) );

    xi = linspace(xi_i, xi_o, M); dxi = abs( xi(2) - xi(1) );

    [chi, wts] = lgwt(Ng, -1,1);
    chi = chi'; wts = wts';

    LegPols = zeros(N_modes+2, Ng);

    for k = 0:N_modes+1

        % NOTE: LegPols(n+1,:) = L_{n}(chi)
        Lk_all = legendre(k, chi);

        LegPols(k+1,:) = Lk_all(1,:);

    end

    I_tensor = I_lin_BiSp_proj(N_modes, xi, chi, wts, LegPols);

    LHS_mat = zeros(N_modes*M); RHS_vec = zeros(N_modes*M,1);

    n_vals = 1:N_modes;

    fn_d2Un = n_vals.^2 + n_vals + 5/4;
    fn_Un   = ( n_vals.^4/2 + n_vals.^3 - n_vals.^2/4 - 3*n_vals/4 + 9/32 );

    coeffs_Un_xi_pm_2dxi = -1/2/dxi^4*ones(1,N_modes);
    coeffs_Un_xi_pm_1dxi = 2/dxi^4 + fn_d2Un/dxi^2;

    coeffs_Un_xi = -(3/dxi^4 + 2/dxi^2*fn_d2Un + fn_Un);

    for l = 3:M-2

        for i = 1:N_modes

            u_xi_l_inds = l:M:(N_modes-1)*M+l;

            %%% number of rows of nZeroCols is always 5
            %%% number of columns of nZeroCols is N_modes

            nZeroCols = [u_xi_l_inds-2; ...
                u_xi_l_inds-1; ...
                u_xi_l_inds; ...
                u_xi_l_inds+1; ...
                u_xi_l_inds+2];

            %%% each row of nZeroCols contains col. no. of U_n(xi_l) with 1 <= n <= N_modes

            %%% each column of nZeroCols contains col. no. of U_n(xi_l) at
            %%% xi values going from xi_l-2*dxi to xi_l+2*dxi

            nZeroRow = (i-1)*M + l; % this is the g.d.e. for U_i(xi_l) for which the current row is being populated

            I_lin = I_tensor(l,i,1:N_modes); I_lin = I_lin(:); I_lin = I_lin';

            LHS_mat( nZeroRow, nZeroCols(1,:) ) = I_lin.*coeffs_Un_xi_pm_2dxi;
            LHS_mat( nZeroRow, nZeroCols(5,:) ) = I_lin.*coeffs_Un_xi_pm_2dxi;

            LHS_mat( nZeroRow, nZeroCols(2,:) ) = I_lin.*coeffs_Un_xi_pm_1dxi;
            LHS_mat( nZeroRow, nZeroCols(4,:) ) = I_lin.*coeffs_Un_xi_pm_1dxi;

            LHS_mat( nZeroRow, nZeroCols(3,:) ) = I_lin.*coeffs_Un_xi;

            [th_i, ~] = theta_modes(i, abs(xi(l)), abs(xi_i), abs(xi_o));

            [th_im1, dth_im1] = theta_modes(i-1, abs(xi(l)), abs(xi_i), abs(xi_o));
            [th_ip1, dth_ip1] = theta_modes(i+1, abs(xi(l)), abs(xi_i), abs(xi_o));

            if i ~= N_modes
                RHS_vec( nZeroRow ) = 2*i*(i+1)/(2*i+1)*( th_i - ...
                    cosh( xi(l) )/2*( th_im1 + th_ip1 ) - ...
                    sinh( abs(xi(l)) )*( dth_im1/(2*i-1) - dth_ip1/(2*i+3) ) );
            elseif i == N_modes
                RHS_vec( nZeroRow ) = 2*i*(i+1)/(2*i+1)*( th_i - ...
                    cosh( xi(l) )/2*( th_im1 ) - ...
                    sinh( abs(xi(l)) )*( dth_im1/(2*i-1) ) );
            end

        end

    end

    for l = [1, 2, M-1, M]

        for n = 1:N_modes

            rowNo = (n-1)*M + l;

            if l == 1 || l == M

                LHS_mat( (n-1)*M+l, (n-1)*M+l ) = 1;

            elseif l == 2

                % minus of the usual because xi has been defined from large to small values
                LHS_mat( rowNo, (n-1)*M + 1 ) = 3/2;
                LHS_mat( rowNo, (n-1)*M + 2 ) = -2;
                LHS_mat( rowNo, (n-1)*M + 3 ) = 1/2;

            else

                % minus of the usual because xi has been defined from large to small values
                LHS_mat( rowNo, (n-1)*M + M     ) = -3/2;
                LHS_mat( rowNo, (n-1)*M + M - 1 ) = 2;
                LHS_mat( rowNo, (n-1)*M + M - 2 ) = -1/2;

            end

        end

    end

    LHS_mat = LHS_mat*dxi^4;
    RHS_vec = RHS_vec*a^4/2*dxi^4;

    Un = LHS_mat\RHS_vec; % U_n(xi) goes from (n-1)*M+1 to n*M

    [~, chi_plot_vals] = chi_vals_surf(e,kappa, N_chi_plot, 'o');

    chi_plot_vals(1) = 1-1e-8; chi_plot_vals(end) = -1+1e-8;
    if e < 0
        chi_plot_vals(1) = 1-1e-8; chi_plot_vals(end) = -1+1e-8;
    end
    
    [uz_max, uz_min, ~] = plot_results(Un, xi, N_modes, chi_plot_vals, kappa, e, plot_vel);

else

    fprintf('\nError: nuclear radius + eccentricity = %.2f >= 1; not possible!\n\n', kappa + abs(e))

end

end