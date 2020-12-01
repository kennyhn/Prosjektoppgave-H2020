clc; close all;
%% Create all constants
constants

%% Disturbance
% Current (disturbance). Constant
%V_x     = 0.25; % m/s in x-axis inertial frame
%V_y     = -0.6; % m/s in y-axis inertial frame
V_x_mat = load('Vx_disturbance.mat').u;
V_y_mat = load('Vy_disturbance.mat').v;
g_z     = 0.91; % Restoring forces. Slightly buoyant

%% References
psi_r1  = -45; %deg
psi_r2  = -45; %deg
u_r     = 0.2; %m/s
v_r     = 0; %m/s
z_r     = 1; %m

% Guidance law parameters
Delta   = 25; % Lookahead distance
x_start = -2*50;
y_start = 1*50;
x_los   = 2*50;
y_los   = 1*50;

alpha_los = atan2(y_los-y_start,x_los-x_start);

zeta_ref    = 1; % critical damping
omega_ref   = 0.4; % Desired bandwidth
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
t_sim = 668; %s

%% Run simulation
sim_output = sim('simulering_ROV_STA_controller.slx');

%% Parse out results
nu              = sim_output.nu.signals.values;
eta             = sim_output.eta.signals.values;

tau_unsat       = sim_output.tau_unsat.signals.values;
tau_sat         = sim_output.tau_sat.signals.values;

psi_d           = sim_output.psi_d.signals.values;
u_d             = sim_output.u_d.signals.values;
v_d             = sim_output.v_d.signals.values;

Vc              = sim_output.disturbance.signals.values;


Vx              = Vc(:, 1);
Vy              = Vc(:, 2);

% Tilstander
x               = eta(:, 1);
y               = eta(:, 2);
z               = eta(:, 3);
psi             = eta(:, 4);
u               = nu(:, 1);
v               = nu(:, 2);
w               = nu(:, 3);
r               = nu(:, 4);


% Control inputs
tau_u_unsat     = tau_unsat(:, 1);
tau_v_unsat     = tau_unsat(:, 2);
tau_w_unsat     = tau_unsat(:, 3);
tau_r_unsat     = tau_unsat(:, 4);

tau_u_sat       = tau_sat(:, 1);
tau_v_sat       = tau_sat(:, 2);
tau_w_sat       = tau_sat(:, 3);
tau_r_sat       = tau_sat(:, 4);

% Parameterestimater for forstyrrelsene

v_xu_hat        = sim_output.V_xu_hat.signals.values;
v_yu_hat        = sim_output.V_yu_hat.signals.values;

v_xv_hat        = sim_output.V_xv_hat.signals.values;
v_yv_hat        = sim_output.V_yv_hat.signals.values;

v_xr_hat        = sim_output.V_xr_hat.signals.values;
v_yr_hat        = sim_output.V_yr_hat.signals.values;
v_xr_hat_sq     = sim_output.V_xr_hat_sq.signals.values;
v_yr_hat_sq     = sim_output.V_yr_hat_sq.signals.values;
v_xyr_hat       = sim_output.V_xyr_hat.signals.values;

time            = sim_output.eta.time;

%% Run plot script
plot_simulation_feedback