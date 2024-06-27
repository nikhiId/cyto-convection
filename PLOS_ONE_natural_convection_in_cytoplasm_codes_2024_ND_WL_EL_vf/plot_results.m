% this function plots the velocity field and vectors once it accepts following arguments:
% the bi-spherical velocity modes: Un
% the discreized values of the xi-coordinate over which Un has been computed: xi
% the total number of bi-spherical modes: N_modes
% the bi-spherical coordinate chi on which to evaluate the velocity: chi
% the nucleus radius (kappa) and eccentricity (e)
% the velocity component (x or z) that the user wishes to plot: plot_vel
% the user can also plot the magnitude of the velocity instead of one component

% the function also returns the following:
% 1. the global maximum (uz_max) and minimum (uz_min) z-velocities
% 2. the values of the Stokes streamfunction: Psi (eqn. (34))

function [uz_max, uz_min, Psi] = plot_results(Un, xi, N_modes, chi, kappa, e, plot_vel)

M = numel(xi);

xi_i = xi(1);
xi_o = xi(end);

a = abs( sinh(xi_o) );

dxi = abs( xi(1) - xi(2) );

Psi = zeros(numel(xi), numel(chi)); uxi = Psi; uchi = uxi;

ur = uxi;
uz = uxi;

r = Psi; z = r;

for l = 1:numel(xi)

    Gamma = cosh( xi(l) ) - chi;

    G_term = Gamma.^(-3/2);

    for n = 1:N_modes

        if l == 1
            dU_n = 1/dxi*( -3/2*Un( (n-1)*M + 1 ) + 2*Un( (n-1)*M + 2 ) - 1/2*Un( (n-1)*M + 3 ) );
        elseif l == numel(xi)
            dU_n = 1/dxi*(  3/2*Un( n*M ) - 2*Un( n*M - 1 ) + 1/2*Un( n*M - 2 ) );
        else
            dU_n = ( Un( (n-1)*M + l + 1 ) - Un( (n-1)*M + l - 1 ) )/2/dxi;
        end

        if xi_i > xi_o
            dU_n = -dU_n;
        end

        U_n = Un( (n-1)*M + l );

        Ln_all = legendre(n, chi); Lnm1_all = legendre(n-1, chi);

        Ln = Ln_all(1,:); Lnm1 = Lnm1_all(1,:);

        Ln_term_Psi = n*( Lnm1 - chi.*Ln );

        Ln_term_uchi = Ln_term_Psi./sqrt(1 - chi.^2);

        Psi(l,:) = Psi(l,:) + G_term.*Ln_term_Psi*U_n;

        uxi(l,:) = uxi(l,:) + (-1)*( -sqrt(Gamma)/a^2*n*(n+1).*Ln*U_n + 3/2/a^2*U_n*Ln_term_Psi./sqrt(Gamma) );

        uchi(l,:) = uchi(l,:) + (-1)*( -sqrt(Gamma)/a^2.*Ln_term_uchi*dU_n + 3/2/a^2*sinh( xi(l) )*U_n*Ln_term_uchi./sqrt(Gamma) );

        exi_er = -sqrt(1 - chi.^2)*sinh( xi(l) )./Gamma;
        echi_er = ( 1 - chi*cosh( xi(l) ) )./Gamma;

        exi_ez = echi_er;
        echi_ez = -exi_er;

        ur(l,:) = uxi(l,:).*exi_er + uchi(l,:).*echi_er;

        uz(l,:) = uxi(l,:).*exi_ez + uchi(l,:).*echi_ez;

    end

    r(l,:) = a*sqrt(1-chi.^2)./Gamma;
    z(l,:) = a*sinh( xi(l) )./Gamma;

    if l == 1 || l == numel(xi)

        alpha = linspace(0,2*pi,100);

        if l == numel(xi)
            x_sp = kappa*sin(alpha);
            z_sp = a*coth(xi_i) + kappa*cos(alpha);
        else
            figure(), hold on, hold all, box on
            
            x_sp = sin(alpha);
            z_sp = a*coth(xi_o) + cos(alpha);
        end
        plot( x_sp, z_sp, '-', 'color', 'k', 'linewidth', 3 ), hold all
    end

end

u_mag = sqrt(ur.^2 + uz.^2);

qgap_xi = 10; qgap_chi = 5;

r_plot = r(1:qgap_xi:end, 1:qgap_chi:end);
z_plot = z(1:qgap_xi:end, 1:qgap_chi:end);

ur_plot = ur(1:qgap_xi:end, 1:qgap_chi:end);
uz_plot = uz(1:qgap_xi:end, 1:qgap_chi:end);

if strcmp(plot_vel, 'z')
    [~, h] = contourf(r,z, uz);
elseif strcmp(plot_vel, 'x')
    [~, h] = contourf(r,z, ur);
elseif strcmp(plot_vel, 'mag')
    [~, h] = contourf(r,z, u_mag);
else
    fprintf('\nError: velocity component to plot not correctly specified...\n')
end

if exist('h', 'var')
    set(h,'LineColor','none')
    h.LevelStep = 5e-5;
end

if strcmp(plot_vel, 'z')
    [~, h] = contourf(-r,z, uz);
elseif strcmp(plot_vel, 'x')
    [~, h] = contourf(-r,z, -ur);
elseif strcmp(plot_vel, 'mag')
    [~, h] = contourf(-r,z, u_mag);
else
    fprintf('\nError: velocity component to plot not correctly specified...\n')
end

if exist('h', 'var')
    set(h,'LineColor','none')
    h.LevelStep = 5e-5;
end

ax = gca;

if ~strcmp(plot_vel, 'mag')
    ax.CLim = [-0.003 0.003];
else
    ax.CLim = [0 0.0035];
end

quiver(r_plot,z_plot, ur_plot,uz_plot, 0.75, 'color', [0 0 0], 'linewidth', 1 )
quiver(-r_plot,z_plot, -ur_plot,uz_plot, 0.75, 'color', [0 0 0], 'linewidth', 1 )

axis tight
axis equal

colorbar
colormap(jet)

xlabel( '$x$', 'interpreter', 'latex', 'fontsize', 22 );
ylabel( '$z$', 'interpreter', 'latex', 'fontsize', 22 );
set( gca, 'fontsize', 22, 'fontname', 'times new roman' );

if xi_i < xi_o
    titleStr = strcat( '$\kappa=$', num2str(kappa), '$,\;e=$', num2str(e) );
else
    titleStr = strcat( '$\kappa=$', num2str(kappa), '$,\;e=-$', num2str(abs(e)) );
end
title( titleStr, 'Interpreter', 'latex', 'fontsize', 22 )

uz_max = max( uz(:) );
uz_min = min( uz(:) );

end