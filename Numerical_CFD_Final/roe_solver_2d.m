function F_hat_ij = roe_solver_2d(W_i, W_j, nu_ij)

% Normalising Edge Area

    mag = sqrt((nu_ij(1)^2 + nu_ij(2)^2));
    nu_ij_norm(1,:) = nu_ij(1)/mag;
    nu_ij_norm(2,:) = nu_ij(2)/mag;

gam = 1.4;

% Left variables
rho_i = W_i(1); % kg/m^3

vx_i = (W_i(2)./W_i(1)); % m/s
vy_i = (W_i(3)./W_i(1)); % m/s

v_mag_i = sqrt(vx_i^2 + vy_i^2);

p_i = (gam-1)*(W_i(4)-rho_i.*(v_mag_i.^2)/2); % N/m^2
h_i = (((p_i/(gam-1)) + ((rho_i)*(v_mag_i^2)/2)) + p_i)/rho_i;

% Right variables
rho_j = W_j(1); % kg/m^3
vx_j = (W_j(2)./W_j(1)); % m/s
vy_j = (W_j(3)./W_j(1)); % m/s

v_mag_j = sqrt(vx_j^2 + vy_j^2);

p_j = (gam-1)*(W_j(4)-rho_j*(v_mag_j.^2)/2); % N/m^2
h_j = (((p_j/(gam-1)) + ((rho_j)*(v_mag_j^2)/2)) + p_j)/rho_j;

% Compress W_i and W_j

v_i = [vx_i ; vy_i];
vn_i = dot(v_i, nu_ij_norm);

E_comp_i = (1/gam)*rho_i*h_i + (1/2*gam)*(gam - 1)*rho_i*(vn_i)^2;
W_comp_i = [rho_i; rho_i*vn_i; E_comp_i] ;

v_j = [vx_j ; vy_j];
vn_j = dot(v_j, nu_ij_norm);

E_comp_j = (1/gam)*rho_j*h_j + (1/2*gam)*(gam - 1)*rho_j*(vn_j)^2;
W_comp_j = [rho_j; rho_j*vn_j; E_comp_j] ;

v_ij = (((rho_j)^0.5)*vn_j + ((rho_i)^0.5)*vn_i)/(((rho_j)^0.5) + ((rho_i)^0.5));
h_ij = (((rho_j)^0.5)*h_j + ((rho_i)^0.5)*h_i)/(((rho_j)^0.5) + ((rho_i)^0.5));

A_ij = [0, 1, 0; ((gam-3)/2)*(v_ij^2), (3-gam)*v_ij, (gam - 1); -v_ij*h_ij + 0.5*(gam - 1)*(v_ij^3), h_ij - (gam - 1)*(v_ij^2), gam*v_ij];

[W_comp_M] = roesolver2D(W_comp_i, W_comp_j);

rho_c = W_comp_M(1);
vn_com = (W_comp_M(2)/rho_c);

e_riem = (W_comp_M(3) - 0.5*rho_c*(vn_com^2))/rho_c;

v_exp = vn_com*(nu_ij_norm) + 0.5 * ((v_i - vn_i*nu_ij_norm) + (v_j - vn_j*nu_ij_norm));

v_riem_exp_mag = sqrt(v_exp(1)^2 + v_exp(2)^2);

W_exp_M = [rho_c; rho_c*v_exp(1); rho_c*v_exp(2);...
    (rho_c*e_riem + 0.5*rho_c*v_riem_exp_mag^2)];

rho_exp = W_exp_M(1);
p_exp_M = (gam-1)*(W_exp_M(4)-rho_exp*(v_riem_exp_mag^2)/2);

V_exp_M = [rho_exp; v_exp(1); v_exp(2); p_exp_M];

F = [rho_exp*transpose(v_exp); rho_exp*v_exp(1)*transpose(v_exp) + p_exp_M*[1,0];...
 rho_exp*v_exp(2)*transpose(v_exp) + p_exp_M*[0,1]; (W_exp_M(4)+ p_exp_M)*transpose(v_exp)];

F_hat_ij = F*nu_ij;

end