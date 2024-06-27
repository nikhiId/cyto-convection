% values of the bi-spherical coordinate chi, evaluated along uniformly
% placed points on the surface of the inner or outer sphere (determined by
% the variable 'in_or_out')

function [alpha, chi] = chi_vals_surf(e,kappa, N_th, in_or_out)

alpha = linspace(1e-8, pi-1e-8, N_th);

xi_i = -e/abs(e)*acosh( (1 - kappa^2 - abs(e)^2)/2/kappa/abs(e) );
xi_o = -e/abs(e)*acosh( (1 - kappa^2 + abs(e)^2)/2/abs(e) );

a = sinh(xi_o);

if strcmp(in_or_out, 'o')
    r = sin(alpha);
    z = cosh(xi_o) + cos(alpha);
else
    r = kappa*sin(alpha);
    z = sinh(xi_o)*coth(xi_i) + kappa*cos(alpha);
end

R = sqrt(r.^2 + z.^2);

Q = sqrt( (R.^2 + a^2).^2 - (2*a*z).^2 );

chi = (R.^2 - a^2)./Q;

end