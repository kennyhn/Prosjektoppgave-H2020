clc; close all; clear;
%% Set scenario
save_simulation     = 0; % 1 for true 0 for false
step_response       = 1; % 1 for step response 0 for guidance
nonlinear_damping   = 1; % 1 to turn on 0 to turn off
coriolis_effect     = 1; % 1 to turn on 0 to turn off

%% Create all constants
constants

filename            = 'simulation_output/FL_controller/FL_controller';

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
% Velocity controller gains
pole_u      = 6.73611111; % place pole at -value. Tune here
k_p_u       = pole_u-d_11/m_11;
gamma_u     = 2;
pole_v      = 6.234567901; % place pole at -value. Tune here
k_p_v       = pole_v-d_22/m_22;
gamma_v     = 2;


% Depth controller gains
zeta_d_heave = 0.5; % Critical damping
wb_d_heave   = 0.35; % desired bandwidth on heave
wn_heave     = wb_d_heave/sqrt(1-2*zeta_d_heave^2+sqrt(4*zeta_d_heave^4-4*zeta_d_heave^2+2));

k_p_z        = wn_heave^2*m_33; 
k_d_z        = 2*m_33*zeta_d_heave*wn_heave-d_33;
k_i_z        = wn_heave/50*k_p_z;

% Heading controller gains
lambda       = 0.79;
k_p_r        = 0.79;
k_p_psi      = 1.89;

gamma_r      = 2;


%% Simulation parameters
if step_response == 1
    t_sim = 1122; %s
else
    t_sim = 1010; %s
end
%% Run simulation
sim_output = sim('simulering_ROV_FL_controller.slx');


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