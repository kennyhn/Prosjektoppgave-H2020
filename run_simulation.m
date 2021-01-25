clc; close all; clear;

%% Set scenario
save_simulation     = 0; % 1 for true 0 for false
use_saved_file      = 0; % 1 for true 0 for false
step_response       = 0; % 1 for step response 0 for guidance
nonlinear_damping   = 1; % 1 to turn on 0 to turn off
coriolis_effect     = 1; % 1 to turn on 0 to turn off
controller          = 0; % 0 for DP model, 1 for FL controller, 2 for PID, controller, 3 for SM controller, 4 for STA.
%% Create all constants
constants

%% Simulate and plot figures

if use_saved_file == 1
    filename   = 'simulation_output/';
   
   
    if isequal(controller, 0)
        filename = strcat(filename, 'dp_model/dp_model');
    elseif isequal(controller, 1)
        filename = strcat(filename, 'FL_controller/FL_controller');
    elseif isequal(controller, 2)
        filename = strcat(filename, 'PID_controller/PID_controller');
    elseif isequal(controller, 3)
        filename = strcat(filename, 'SM_controller/SM_controller');
    elseif isequal(controller, 4)
        filename = strcat(filename, 'STA_controller/STA_controller');
    end
   
    if step_response == 0
        filename = strcat(filename, '_guidance.mat');
    else  
        filename = strcat(filename, '_step.mat');
    end
   
    sim_output = load(filename).sim_output;
    if isequal(controller, 0)
        plot_simulation_dp
    else
        plot_simulation_feedback
    end
else
    % Simulation parameters
    if step_response == 1
        t_sim = 1122; %s
    else
        t_sim = 1450; %s
    end

    % References
    if step_response == 1
        u_r     = 0.2; % m/s
        v_r     = 0; % m/s
    else
        u_r     = 0;%0.14; %m/s
        v_r     = 0.2;%0.14; %m/s
    end
    psi_r   = deg2rad(-90);%-45);
    z_r     = 10; % m

    psi_r1      = 0; % deg
    psi_r2      = 45; % deg
    time_step   = 700; % seconds. about right after passing the farm
    psi_r1      = deg2rad(psi_r1); % rad
    psi_r2      = deg2rad(psi_r2); % rad


    % Guidance law parameters
    Delta   = 25; % Lookahead distance

    % Reference model parameters
    zeta_ref    = 1; % critical damping
    omega_ref   = 1.5; % Desired bandwidth
    T_ref       = 0.2; % Desired time constant for first-order model

    if isequal(controller, 0) % DP controller
        filename            = 'simulation_output/dp_model/dp_model';
       
        % Run simulation and plotting files
        tuning_dp_model
        sim_output = sim('simulering_ROV_DP_model.slx');
        plot_simulation_dp
    elseif isequal(controller, 1) % FL controller
        filename            = 'simulation_output/FL_controller/FL_controller';
        
        % Run simulation and plotting files
        tuning_FL_controller
        sim_output = sim('simulering_ROV_FL_controller.slx');
        plot_simulation_feedback
    elseif isequal(controller, 2) % PID controller
        filename            = 'simulation_output/PID_controller/PID_controller';
        
        % Run simulation and plotting files
        tuning_PID_controller
        sim_output = sim('simulering_ROV_PID_controller.slx');
        plot_simulation_feedback
    elseif isequal(controller, 3) % SM controller
        filename            = 'simulation_output/SM_controller/SM_controller';
        
        % Run simulation and plotting files
        tuning_SM_controller
        sim_output = sim('simulering_ROV_SM_controller.slx');
        plot_simulation_feedback
    elseif isequal(controller, 4) % STA controller
        filename            = 'simulation_output/STA_controller/STA_controller';
        
        % Run simulation and plotting files
        tuning_STA_controller
        sim_output = sim('simulering_ROV_STA_controller.slx');
        plot_simulation_feedback
    end


    if save_simulation == 1
        if step_response == 1
            filename = strcat(filename, '_step.mat');
        else
            filename = strcat(filename, '_guidance.mat');
        end
        save(filename, 'sim_output');
    end
end
