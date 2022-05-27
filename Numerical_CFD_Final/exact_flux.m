% Exact flux
% Primitive state variable vector V_ij
% normal ν_ij

function F = exact_flux(V_ij, nu_ij)

gam = 1.4;
rho = V_ij(1);
u = V_ij(2);
v = V_ij(2);
p = V_ij(2);
E = 0.5 * rho * (u^2 + v^2) + p / (gam - 1);
v_n = u * nu_ij(1) + v * nu_ij(2);

F = [rho * v_n; rho * u * v_n + p * nu_ij(1); rho * v * v_n + p * nu_ij(2); (E + p) * v_n];
end