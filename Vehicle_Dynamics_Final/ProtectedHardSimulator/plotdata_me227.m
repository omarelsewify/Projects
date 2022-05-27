function plotdata_me227(varargin)
--------------------------------------------------------------------------
% HEADER
--------------------------------------------------------------------------
ME227 data visualization
Modified by John Talbot 4/15/21
Modified by Alaisha Alexander 5/10/21
Modified by Will Harvey 5/12/21
  - Updated subplot(2,2,1)
  - Added a final print statement to indicate code completion

close all; clc;


--------------------------------------------------------------------------
% FORMATTING
--------------------------------------------------------------------------
Set all interpreters to Latex
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');


--------------------------------------------------------------------------
% OPEN DATA
--------------------------------------------------------------------------
if isempty(varargin)
    [LoadFileName, LoadPathName] = uigetfile('*.mat', 'Select a .mat file');
else    
    LoadPathName = '';
    LoadFileName = varargin{1};
end

try    
    load([LoadPathName LoadFileName]);
catch
    disp('No data file was successfully loaded!');
    return;
end

--------------------------------------------------------------------------
% PLOT DATA
--------------------------------------------------------------------------
f = figure('Name', 'ME227 Project Experiment', 'IntegerHandle', 'off');

    Longitudinal speed
    subplot(2,2,1); hold on;
    plot(t_, Ux_ .* ones(1,length(t_)) , 'b--')
    plot(t_, Ux_, 'r-')
    grid on
    xlabel('Time [s]')
    ylabel('Longitudinal Velocity, $u_x$ [m/s]')
    legend('Desired Speed', 'Actual Speed', 'Location', 'South')
    
    Steer Angle
    subplot(2,2,2)
    plot(t_, rad2deg(delta_), 'r-')
    grid on
    xlabel('Time [s]')
    
    ylabel('Steer Angle, $\delta$ [deg]')
    
    Lateral error
    subplot(2,2,3)
    plot(t_, e_, 'r-')
    grid on
    xlabel('Time [s]')
    ylabel('Lateral Error, $e$ [m]')
    
    Longitudinal Force
    subplot(2,2,4)
    plot(t_, Fx_, 'r-')
    grid on
    xlabel('Time [s]')
    ylabel('Longitudinal Force, $F_x$ [N]')

disp("Script complete, but the plot may not appear for another minute or two.");