%--------------------------------------------------------------------------
%% HEADER
%--------------------------------------------------------------------------
% ME227 Spr 2021
% Homework 4 - Question 2.B
close all; clear; clc;
load('project_path.mat')
%--------------------------------------------------------------------------
%% CONSTANTS AND PARAMS
%--------------------------------------------------------------------------
% Load Niki params
setup_niki;
% Gains and conditions
K_la = 3000;        % [N/m]
x_la_ = 0:2:20;     % [m]    
Ux = mean(path.UxDes);            % [m/s]
lenxla = length(x_la_);

% Allocate space for poles (We know there are 4)
poles_ = zeros(4,lenxla);

%--------------------------------------------------------------------------
%% CREATE SYSTEM MATRIX
%--------------------------------------------------------------------------
for idx = 1:lenxla
    % Select speed
    x_la = x_la_(idx);
    
    %{
    To make the state matrix easier to input, create each term separately
    here according to this template - we'll complile these into the matrix
    at the end. We recommend you keep this general and let MATLAB fill in
    each of the values as you set them up above. Then you can copy and
    paste this section into later problems.
    
        A = [aM,  bM,  cM,  dM]
            [eM,  fM,  gM,  hM]
            [iM,  jM,  kM,  lM]
            [mM,  nM,  oM,  pM]
    %}
    
    aM = 0;
    bM = 1;
    cM = 0;
    dM = 0;
    eM = -K_la/veh.m;
    fM = -(f_tire.Ca_lin+r_tire.Ca_lin)/(veh.m*Ux);
    gM = ((f_tire.Ca_lin+r_tire.Ca_lin)/(veh.m))-K_la*x_la/veh.m;
    hM = (-veh.a*f_tire.Ca_lin+veh.b*r_tire.Ca_lin)/(veh.m*Ux);
    iM = 0;
    jM = 0;
    kM = 0;
    lM = 1;
    mM = -K_la*veh.a/veh.Iz;
    nM = (veh.b*r_tire.Ca_lin-veh.a*f_tire.Ca_lin)/(veh.Iz*Ux);
    oM = ((veh.a*f_tire.Ca_lin-veh.b*r_tire.Ca_lin)/(veh.Iz))-(K_la*x_la*veh.a)/(veh.Iz);
    pM = -((veh.a^2)*f_tire.Ca_lin+(veh.b^2)*r_tire.Ca_lin)/(veh.Iz*Ux);

    A = [[aM,  bM,  cM,  dM];
         [eM,  fM,  gM,  hM];
         [iM,  jM,  kM,  lM];
         [mM,  nM,  oM,  pM]];
    
   % Calculate pole positions
   poles_(:,idx) = eig(A);
end

%--------------------------------------------------------------------------
%% PLOT RESULTS
%--------------------------------------------------------------------------
figure
cmap = colormap(winter(lenxla));
for idx = 1:lenxla
    plot(real(poles_(:,idx)), imag(poles_(:,idx)), 'x', 'Color', cmap(idx,:))
    hold on
end
xline(0,'--');
grid on
xlabel('Real Axis')
ylabel('Imaginary Axis')
cbar = colorbar('Ticks', x_la_);
caxis([x_la_(1) x_la_(end)])
cbar.Label.String = 'x_{la} [m]';

