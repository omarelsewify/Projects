function [ delta, Fx ] = groupnum_controller( s, e, dpsi, Ux, Uy, r, control_mode, path)
%ME227 Controller:
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


%--------------------------------------------------------------------------
%% Control Parameters
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
%% Lateral Control Law
%--------------------------------------------------------------------------
%Use the Lateral Control Law to Calculate Delta

if control_mode == 1 %Lookahead Controler
    
else %Your second controller

end

%--------------------------------------------------------------------------
%% Longitudinal Control Law
%--------------------------------------------------------------------------
%Use the Longitudinal Control Law to Calcuate Fx


end
