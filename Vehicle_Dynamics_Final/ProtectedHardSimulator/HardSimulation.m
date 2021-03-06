%--------------------------------------------------------------------------
%% HEADER
%--------------------------------------------------------------------------
% ME227 Spr 2021
% Nonlinear Vehicle Hard Simulation - Project

clear; clc; close all;

table_path = "/Users/omarelsewify/Google Drive/MS MechEng/AA222 Engineering Design Optimisation/Final_Project/Sim3_Results.xlsx";

PID_results = readmatrix(table_path);
PID_results = PID_results(PID_results(:,1) == 3,:);

for PID_row_idx = 1:length(PID_results)

    PID_row = PID_results(PID_row_idx,:);
    
    % Only PID gains with no existing output
    if ~isnan(PID_row(5))
        continue
    end

    simtype = PID_row(1);
    K_la = PID_row(2);
    x_la = PID_row(3);
    K_long = PID_row(4);

    % Clear persistent variables
    clear simulate_step; clear gps_model;
    clear brakeDynamics; clear engineDynamics; clear steeringDynamics;


    %--------------------------------------------------------------------------
    %% CONTROLLER SETUP
    %--------------------------------------------------------------------------
    % Select lookahead controller mode
    % control_mode = 1; % Lookahead Controller
    control_mode = 2; % Your second controller!

    %--------------------------------------------------------------------------
    %% CONSTANTS AND PARAMS
    %--------------------------------------------------------------------------
    % Load Niki params
    setup_niki;

    % Create time vector
    dt = 0.001;
    t_ = 0:dt:40; %Note: 35 seconds is about how long it takes to traverse the path
    lenT = length(t_);

    % Load path and speed profile
    load('project_path.mat')

    % Straight Path
    % path.k_1pm = zeros(length(path.k_1pm),1);
    % path.posE_m = zeros(length(path.s_m),1);
    % path.psi_rad = zeros(length(path.k_1pm),1);

    % Set simulation mode
    % sim_mode = 0;   % Nonlinear model from homework
    % sim_mode = 1;   % Actuator dynamics only
    % sim_mode = 2;   % Actuator dynamics and noise on measurements
    sim_mode = 3;   % Actuator dynamics, noisy measurements, and hold period

    wait1 = waitbar(0, 'Simulation Initializing');
    j = 0;

    %--------------------------------------------------------------------------
    %% ALLOCATE MEMORY
    %--------------------------------------------------------------------------
    % Allocate space for results (!!Do not change these variable names - store
    % your results in these variables!!)
    s_ = zeros(1,lenT);
    e_ = zeros(1,lenT);
    dpsi_ = zeros(1,lenT);
    r_ = zeros(1,lenT);
    Ux_ = zeros(1,lenT);
    Uy_ = zeros(1,lenT);
    delta_ = zeros(1,lenT);
    Fx_ = zeros(1,lenT);
    Fy_f = zeros(1,lenT);
    Fy_r = zeros(1,lenT);
    Fx_f = zeros(1,lenT);
    Fx_r = zeros(1,lenT);
    Fyf_max = zeros(1,lenT);
    Fyr_max = zeros(1,lenT);
    ax_ = zeros(1,lenT);
    ay_ = zeros(1,lenT);
    atot_ = zeros(1,lenT);
    delta_actual_ = zeros(1,lenT);

    % Allocate space for UxDes_
    UxDes_ = zeros(1, lenT);

    %--------------------------------------------------------------------------
    %% SET INITIAL CONDITIONS
    %--------------------------------------------------------------------------
    s_(1) = 0.5;
    e_(1) = 0.15;
    dpsi_(1) = 0;
    r_(1) = 0;
    Ux_(1) = 1;
    Uy_(1) = 0;

    %--------------------------------------------------------------------------
    %% SIMULATION LOOP
    %--------------------------------------------------------------------------
    k = 0;
    % Loop through every time step
    for i = 1:lenT-1
        % look up kappa
        kappa = interp1(path.s_m, path.k_1pm, s_(i));

        % look up UxDes
        UxDes_(i) = interp1(path.s_m, path.UxDes, s_(i));

        if k == 0
            % Calculate actuator commands
            [delta_(i), Fx_(i)] = FINAL_group2_controller(s_(i), e_(i), dpsi_(i), Ux_(i),...
                Uy_(i), r_(i), control_mode, path, ...
                K_la, x_la, K_long);
            k = 4;
        else %First order hold for 4 ms
            delta_(i) = delta_(i-1);
            Fx_(i) = Fx_(i-1);
            k = k - 1;
        end

        if (t_(i) < 5) && (sim_mode == 3)
            % Hold here
            [Ux_(i+1), Uy_(i+1), r_(i+1), s_(i+1), e_(i+1), dpsi_(i+1)] =...
                simulate_hold(delta_(i), Fx_(i), kappa, dt, veh, f_tire, r_tire, 1);
        else
            % Take simulation step
            [Ux_(i+1), Uy_(i+1), r_(i+1), s_(i+1), e_(i+1), dpsi_(i+1),...
                delta_actual_(i), Fx_f(i), Fx_r(i), Fy_f(i), Fy_r(i),...
                Fyf_max(i), Fyr_max(i), ax_(i), ay_(i), atot_(i)] =...
                simulate_step(delta_(i), Fx_(i), kappa, dt, veh, f_tire, r_tire, sim_mode);
        end

        j = j + 1;
        if j > 100
            waitbar(i/(lenT-1), wait1, 'Simulation Running');
            j = 0;
        end

        % If Ux is below zero mps, stop the simulation and plot what we have
        % simulated so far
        if Ux_(i+1) < -0.1
            s_((i+1):end) = s_(i);
            e_((i+1):end) = e_(i);
            dpsi_((i+1):end) = dpsi_(i);
            delta_((i+1):end) = delta_(i);
            delta_actual_((i+1):end) = delta_actual_(i);
            plot_simulation(t_, Ux_, Uy_, r_, s_, e_, dpsi_, delta_, Fx_, Fx_r, Fx_f, ...
                Fy_f, Fyf_max, Fy_r, Fyr_max, ax_, ay_, atot_, ...
                delta_actual_, path, veh)
            disp('The vehicle came to a stop in the simulation.')
            % Display lateral error
            fprintf('Max Lateral Error at CG: %.4f\n', max(e_));
            % Display difference in desired and actual speeds
            fprintf('Max Difference in Desired and Actual Speed: %.4f\n', max(abs(UxDes_(6000:end-1500) - Ux_(6000:end-1500))));
            close(wait1);
            return
        end
    end

    close(wait1);

    %--------------------------------------------------------------------------
    %% PLOT
    %--------------------------------------------------------------------------
%     plot_simulation(t_, Ux_, Uy_, r_, s_, e_, dpsi_, delta_, Fx_, Fx_r, Fx_f, ...
%                 Fy_f, Fyf_max, Fy_r, Fyr_max, ax_, ay_, atot_, ...
%                 delta_actual_, path, veh)

    fprintf('Max Lateral Error at CG: %.4f\n', max(e_));
    % Display difference in desired and actual speeds
    fprintf('Max Difference in Desired and Actual Speed: %.4f\n', max(abs(UxDes_(6000:end-1500) - Ux_(6000:end-1500))));
    % Display max accelerations
    fprintf('Max absolute lateral acceleration ay: %.4f\n', max(abs(ay_)))
    fprintf('Max absolute long acceleration ay: %.4f\n', max(abs(ax_)))
    
    res2add = [3, K_la, x_la, K_long, max(e_), max(abs(UxDes_(6000:end-1500) - Ux_(6000:end-1500))), max(abs(ay_)), max(abs(ax_))];
    add_pos = sprintf('A%i:H%i',PID_row_idx+1,PID_row_idx+1);
    writematrix(res2add, table_path,'Sheet',1,'Range',add_pos)

    uxerror = zeros(size(Ux_));

    for i = 6000:33501
        uxerror(i) = abs(UxDes_(i)-Ux_(i));
    end

end
