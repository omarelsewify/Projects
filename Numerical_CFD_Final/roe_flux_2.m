% Roe's Flux method 2
% Primitive state variable vector V_i, Vj
% normal ν_ij, from node i to node j

function F =roe_flux_2(W_i, W_j, nu_ij)
gam = 1.4;

l = sqrt(nu_ij(1)^2 + nu_ij(2)^2);

% construct unit normal/tagent vector
nu_ij = nu_ij / l;
t_ij = [-nu_ij(2), nu_ij(1)];

% left state
rho_l = W_i(1);
u_l = W_i(2)/W_i(1);
v_l = W_i(3)/W_i(1);
vmagl = sqrt(u_l^2 + v_l^2);
p_l = (gam-1)*(W_i(4)-rho_l*(vmagl^2)/2);
un_l = u_l * nu_ij(1) + v_l * nu_ij(2);
ut_l = u_l * t_ij(1) + v_l * t_ij(2);
w_l = [rho_l, rho_l * un_l, rho_l * (un_l * un_l) / 2 + p_l / (gam - 1)];
h_l = 0.5 * (un_l * un_l) + gam * p_l / (rho_l * (gam - 1.0));

% right state
rho_r = W_j(1);
u_r = W_j(2)/W_j(1);
v_r = W_j(3)/W_j(1);
vmagr = sqrt(u_r^2 + v_r^2);
p_r = (gam-1)*(W_j(4)-rho_r*(vmagr^2)/2);
un_r = u_r * nu_ij(1) + v_r * nu_ij(2);
ut_r = u_r * t_ij(1) + v_r * t_ij(2);
w_r = [rho_r, rho_r * un_r, rho_r * (un_r * un_r) / 2 + p_r / (gam - 1)];
h_r = 0.5 * (un_r * un_r) + gam * p_r / (rho_r * (gam - 1.0));

% compute the Roe-averaged quatities
RT = sqrt(rho_r / rho_l);
rho_rl = RT * rho_l;
un_rl = (un_l + RT * un_r) / (1.0 + RT);
h_rl = (h_l + RT * h_r) / (1.0 + RT);
a_rl = sqrt((gam - 1) * (h_rl - 0.5 * (un_rl * un_rl)));

% wave strengths
du = un_r - un_l;
dp = p_r - p_l;
drho = rho_r - rho_l;
du = [du + dp/(rho_rl * a_rl),  drho - dp/(a_rl * a_rl), du - dp/(rho_rl * a_rl)];

% compute the Roe-average wave speeds
ws = [abs(un_rl + a_rl); abs(un_rl);  abs(un_rl - a_rl)];

%compute the right characteristic vectors
r_1 =  rho_rl/(2 * a_rl) * [1; un_rl + a_rl; h_rl + a_rl * un_rl];
r_2 =  [1.0; un_rl; un_rl * un_rl/2];
r_3 = -rho_rl/(2 * a_rl) * [1; un_rl - a_rl; h_rl - a_rl * un_rl];

%compute the riemannann state
w_riemann = zeros(3);
M = un_rl / a_rl;

if(M <=-1)
    w_riemann = w_r;
elseif(M <= 0 && M >= -1)
    w_riemann = w_r - r_1*du(1);
elseif(M >= 0 && M <= 1)
    w_riemann = w_l + r_3*du(3);
else
    w_riemann = w_l;
end
rho_riemann = w_riemann(1);
un_riemann = w_riemann(2)/w_riemann(1);
p_riemann = (w_riemann(3) - w_riemann(2)^2/w_riemann(1)/2.0) * (gam - 1);

%update vtriemann with upwind scheme;
u_riemann = (un_riemann*nu_ij + (ut_l + ut_r)/2.0*t_ij);
ux_riemann = u_riemann(1);
uy_riemann = u_riemann(2);
V_riemann = [rho_riemann; ux_riemann; uy_riemann; p_riemann];

F = exact_flux(V_riemann, nu_ij*l);
end