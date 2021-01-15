clc; close all;
clc; close all; clear;
%% Set scenario
save_simulation     = 0; % 1 for true 0 for false
step_response       = 1; % 1 for step response 0 for guidance
nonlinear_damping   = 1; % 1 to turn on 0 to turn off
coriolis_effect     = 1; % 1 to turn on 0 to turn off

%% Create all constants
constants

filename            = 'simulation_output/STA_controller/STA_controller';

%% Disturbance
% Current (disturbance). Constant
%V_x     = 0.25; % m/s in x-axis inertial frame
%V_y     = -0.6; % m/s in y-axis inertial frame
V_x_mat = load('Vx_disturbance.mat').u;
V_y_mat = load('Vy_disturbance.mat').v;
g_z     = 0.91; % Restoring forces. Slightly buoyant

%% References
if step_response == 1
    u_r     = 0.2; % m/s
    v_r     = 0; % m/s
else
    u_r     = 0.14; %m/s
    v_r     = 0.14; %m/s
end

psi_r   = deg2rad(-45);
z_r     = 10; %m

psi_r1      = 0; % deg
psi_r2      = 45; % deg
time_step   = 700; % seconds. about right after passing the farm
psi_r1      = deg2rad(psi_r1); % rad
psi_r2      = deg2rad(psi_r2); % rad

% Guidance law parameters
Delta   = 25; % Lookahead distance
x_start = -2*50;
y_start = 1*80;
x_los   = 2*50;
y_los   = 1*80;

alpha_los = atan2(y_los-y_start,x_los-x_start);

zeta_ref    = 1; % critical damping
omega_ref   = 1.5; % Desired bandwidth
T_ref       = 0.2; % Desired time constant for first-order model

%% Controller gains
alpha_m_u   = 1.5;
omega_u     = 0.01;
gam_u       = 1;
epsilon_u   = 0.0001;
lambda_u    = 10;

alpha_m_v   = 1.5;
omega_v     = 0.01;
gam_v       = 1;
epsilon_v   = 0.0001;
lambda_v    = 10;

alpha_m_w   = 1.5;
omega_w     = 0.01;
gam_w       = 1;
epsilon_w   = 0.0001;
lambda_w    = 0.5;

alpha_m_psi = 1.5;
omega_psi   = 0.001;
gam_psi     = 1;
epsilon_psi = 0.0001;
lambda_psi  = 0.5;


%% Adaptive gains
gamma_u     = 0.1;
gamma_v     = 0.1;
gamma_r     = 0.1;
lambda      = 1;

%% Simulation parameters
if step_response == 1
    t_sim = 25; %s
else
    t_sim = 668; %s
end
%% Run simulation
sim_output = sim('simulering_ROV_STA_controller.slx');

if save_simulation == 1
    if step_response == 1
        filename = strcat(filename, '_step.mat');
    else
        filename = strcat(filename, '_guidance.mat');
    end
    save(filename, 'sim_output');
end

%% Run plot script
plot_simulation_feedback