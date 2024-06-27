% evaluation of the I_lin matrix of eqn. (42) in Appendix

function I_matrix = I_lin_BiSp_proj(N_modes, xi, chi, wts, LegPols)

I_matrix = zeros(numel(xi), N_modes, N_modes);

for l = 1:numel(xi)

    Gamma_cube = (cosh( xi(l) ) - chi).^3;

    for i = 1:N_modes

        % NOTE: LegPols(i,:) = L_{i-1}(chi)
        Li_term = (chi.*LegPols(i+1,:) - LegPols(i,:))./( (chi.^2 - 1) ); % .*sqrt(1 - chi.^2)

        for n = 1:N_modes
            
            % NOTE: LegPols(n+1,:) = L_{n}(chi)
            Ln_term = LegPols(n+2,:) - chi.*LegPols(n+1,:);

            integrand = i*(n+1)*Gamma_cube.*Li_term.*Ln_term;

            I_matrix(l,i,n) = sum( integrand.*wts );

        end

    end

end

end