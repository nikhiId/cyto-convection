% given xi_i and xi_o (values of the bi-spherical coordinate xi, for the
% nuclear and cell membrane), this function outputs the n-th temperature
% mode (eqn. (28) in Appendix) and its derivative

function [th_n, dth_n] = theta_modes(n, xi, xi_i, xi_o)

D_n = sinh( (n+1/2)*(xi_i - xi_o) );

a_n = sqrt(2)*exp( -(n+1/2)*xi_i )*cosh( (n+1/2)*xi_o )/D_n;
b_n =-sqrt(2)*exp( -(n+1/2)*xi_i )*sinh( (n+1/2)*xi_o )/D_n;

th_n  = a_n*sinh( (n+1/2)*xi ) + b_n*cosh( (n+1/2)*xi );

dth_n = (n+1/2)*a_n*cosh( (n+1/2)*xi ) + ...
    (n+1/2)*b_n*sinh( (n+1/2)*xi );

end