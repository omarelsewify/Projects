function [ delta, Fx ] = group2_controller( s, e, dpsi, Ux, Uy, r, control_mode, path)
% ME227 Group 2 Controller:
% Spring 2021
% Prof. Chris Gerdes & CAs Nathan Spielberg, John Alsterda, Alaisha
% Alexander, Will Harvey, Lucio Mondavi, John Talbot, Trey Weber
% 

%--------------------------------------------------------------------------
%% Constants
%--------------------------------------------------------------------------
g = 9.81;                       % [m/s^2]  gravity

%--------------------------------------------------------------------------
%% Vehicle Parameters
%--------------------------------------------------------------------------
m  = 1776;                  % [kg]     mass with 2 occupants
Iz = 2763.49;               % [kg-m^2] rotational inertia
a  = 1.264;                 % [m]      distance from CoM to front axle
b  = 1.367;                 % [m]      distance from C0M to rear axle
L  = a + b;         % [m]      wheelbase
Wf = m*g*(b/L); % [N]      static front axle weight
Wr = m*g*(a/L); % [N]      static rear axle weight

%--------------------------------------------------------------------------
%% Tire Parameters
%--------------------------------------------------------------------------
% Front tires
f_tire.Ca_lin = 80000;          % [N/rad]  linear model cornering stiffness
f_tire.Cy     = 110000;         % [N/rad]  fiala model cornering stiffness
f_tire.mu_s   = 0.90;           %          sliding friction coefficient
f_tire.mu     = 0.90;           %          peak friction coefficient

% Rear tires
r_tire.Ca_lin = 120000;
r_tire.Cy     = 180000;
r_tire.mu_s   = 0.94;
r_tire.mu     = 0.94;

%--------------------------------------------------------------------------
%% Find Path Dependent Parameters
%--------------------------------------------------------------------------
% Interpolate desired speed from path
Ux_des = interp1(path.s_m, path.UxDes, s);

% Interpolate desired acceleration from path
ax_des = interp1(path.s_m, path.axDes, s);

% Find Curvature for the current distance along the path via interpolation
kappa = interp1(path.s_m, path.k_1pm, s);

%--------------------------------------------------------------------------
%% Control Parameters
%--------------------------------------------------------------------------
% Gains
K_la = 3000; % 1750
x_la = 8.1; % 19
K_long = 0.158 * g * m;

% Longitudinal parameters
f_rr = 0.015; % rolling resistance friction
C_DA = 0.594; % coefficient of drag + area
rho = 1.225; % air density

% Understeer gradient
K_radpmps2 = (1 / g) * (Wf / f_tire.Ca_lin - Wr / r_tire.Ca_lin); % [rad/m/s^2]

%--------------------------------------------------------------------------
%% Lateral Control Law
%--------------------------------------------------------------------------
% Use the Lateral Control Law to Calculate Delta
if control_mode == 1 % Lookahead Controler
    % Feedback + feedforward
    dpsi_ss = kappa * ((m * a * Ux^2) / (L * r_tire.Ca_lin) - b);
    delta_ff = ((K_la * x_la) / f_tire.Ca_lin) * dpsi_ss + kappa * (L + K_radpmps2 * Ux^2);
    delta = - (K_la / f_tire.Ca_lin) * (e + x_la * dpsi) + delta_ff;
else %Your second controller
    delta = 0;
end

%--------------------------------------------------------------------------
%% Longitudinal Control Law
%--------------------------------------------------------------------------
% Use the Longitudinal Control Law to Calcuate Fx
F_rr = f_rr * m * g; % rolling resistance
F_d = 0.5 * rho * C_DA * Ux^2; % aerodynamic drag
Fx_feedback = K_long * (Ux_des - Ux); % feedback controller
Fx = m * ax_des + F_rr + F_d + Fx_feedback;

end
