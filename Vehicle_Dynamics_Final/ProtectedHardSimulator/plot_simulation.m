function plot_simulation(t_, Ux_, Uy_, r_, s_, e_, dpsi_, delta_, Fx_, Fx_r, Fx_f, ...
            Fy_f, Fyf_max, Fy_r, Fyr_max, ax_, ay_, atot_, ...
            delta_actual_, path, veh)
% Plot all of the vehicle states
figure

% Plot Ux
subplot(2,3,1)
plot(s_, Ux_)
hold on
plot(path.s_m, path.UxDes)
grid on
legend('Simulated', 'Desired')
xlabel('Distance [m]')
ylabel('Ux [m/s]')

% Plot Uy
subplot(2,3,2)
plot(t_, Uy_)
grid on
title('Vehicle States')
xlabel('Time [s]')
ylabel('Uy [m/s]')

% Plot r
subplot(2,3,3)
plot(t_, rad2deg(r_))
grid on
xlabel('Time [s]')
ylabel('r [deg/s]')

% Plot s
subplot(2,3,4)
plot(t_, s_)
grid on
xlabel('Time [s]')
ylabel('s [m]')

% Plot e
subplot(2,3,5)
plot(t_, e_)
grid on
xlabel('Time [s]')
ylabel('e [m]')

% Plot dpsi
subplot(2,3,6)
plot(t_, rad2deg(dpsi_))
grid on
xlabel('Time [s]')
ylabel('dpsi [deg]')


% Plot the actuator commands
figure

subplot(2,1,1)
plot(t_, rad2deg(delta_))
hold on
plot(t_, rad2deg(delta_actual_))
grid on
xlabel('Time [s]')
ylabel('Steer Angle [deg]')
legend('Commanded','Actual')
title('Actuator Commands')

subplot(2,1,2)
plot(t_, Fx_)
hold on
plot(t_, Fx_r + Fx_f)
ylim([-10000,10000]);
grid on
xlabel('Time [s]')
legend('Commanded','Actual')
ylabel('Total Longitudinal Force [Fx]')

% Plot tire forces
figure

subplot(2,1,1)
plot(t_, Fy_f)
hold on
plot(t_, Fyf_max, 'k--')
plot(t_, -Fyf_max, 'k--')
ylim([-1.5e4,1.5e4]);
grid on
xlabel('Time [s]')
ylabel('Front Lateral Force [N]')
title('Front Lateral Force [F_{yf}]')

subplot(2,1,2)
plot(t_, Fy_r)
hold on
plot(t_, Fyr_max, 'k--')
plot(t_, -Fyr_max, 'k--')
ylim([-2e4,2e4]);
grid on
xlabel('Time [s]')
ylabel('Rear Lateral Force [N]')
ylabel('Rear Lateral Force [F_{yr}]')

% Plot Accelerations
figure
plot(s_, ax_)
hold on
plot(s_, ay_)
plot(s_, atot_)
ylim([-8, 8])
grid on
xlabel('Distance Along Path [m]')
ylabel('Acceleration [m/s^2]')
title('Acceleration vs. Path Position')
leg1 = legend('a_x','a_y','a_{tot}');
set(leg1,'Interpreter','tex')

% Plot speed profile
figure
a(1) = subplot(2,1,1); 
plot( path.s_m, path.UxDes, 'linew', 2);
ylabel('U_x, m/s'); 
ylim([-1 15]); grid on
title('Speed Profile')
a(2) = subplot(2,1,2); 
plot(path.s_m, path.axDes, path.s_m, path.ayDes, path.s_m, sqrt(path.axDes.^2 + path.ayDes.^2), 'linew', 2);
linkaxes(a, 'x');  
ylabel('acceleration, m/s^2'); 
xlabel('s, m');
ylim([-4, 5]); 
grid on; 
legend('a_x', 'a_y', 'a_{tot}', 'Location', 'best')
title('Acceleration Profile')

%--------------------------------------------------------------------------
%% ANIMATE VEHICLE
%--------------------------------------------------------------------------
%Note, press 'P' to play, 'M' to record a movie, or mousewheel/arrow keys
%to scroll through simulation
%animate(path, veh, dpsi_, s_, e_, delta_) % You can press p on your keyboard to start the animation

end