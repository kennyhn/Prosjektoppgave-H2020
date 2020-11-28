clc; close all;
%% Create all constants
constants

%% Disturbance
% Current (disturbance). Constant
%V_x     = 0.25; % m/s in x-axis inertial frame
%V_y     = -0.6; % m/s in y-axis inertial frame
V_x_mat = load('Vx_disturbance.mat').u;
V_y_mat = load('Vy_disturbance.mat').v;
g_z     = 0.91; % Restoring forces. Slightly buyoant

%% References
psi_r1  = 0; %deg
psi_r2  = 0; %deg
u_r     = 0.2; %m/s
v_r     = 0; %m/s
z_r     = 10; %m

zeta_ref    = 1; % critical damping
omega_ref   = 0.4; % Desired bandwidth
T_ref       = 0.2; % Desired time constant for first-order model

%% Controller gains
% Velocity controller gains
pole_u      = 5; % place pole at -value. Tune here
k_p_u       = pole_u-d_11/m_11;
pole_v      = 5; % place pole at -value. Tune here
k_p_v       = pole_v-d_22/m_22;

% Depth controller gains
zeta_d_heave = 1; % Critical damping
wb_d_heave   = 0.3; % desired bandwidth on heave
wn_heave     = wb_d_heave/sqrt(1-2*zeta_d_heave^2+sqrt(4*zeta_d_heave^4-4*zeta_d_heave^2+2));

k_p_w        = wn_heave^2*m_33; 
k_d_w        = 2*m_33*zeta_d_heave*wn_heave-d_33;
k_i_w        = wn_heave/10*k_p_w;

% Heading controller gains
zeta_d_psi   = 1; % Critical damping
wb_d_psi     = 0.8; % desired bandwidth on heading
wn_psi       = wb_d_psi/sqrt(1-2*zeta_d_psi^2+sqrt(4*zeta_d_psi^4-4*zeta_d_psi^2+2));

k_p_psi      = wn_psi^2*m_66;
k_d_psi      = 2*m_66*zeta_d_psi*wn_psi-d_66;
k_i_psi      = 1/10*k_p_psi;
%% Adaptive controller gains
gamma1      = 0.5;
gamma2      = 0.5;

%% Simulation parameters
t_sim = 999; %s

%% Run simulation
sim_output = sim('simulering_ROV_DP_model.slx');

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
v_x_hat         = sim_output.V_x_hat.signals.values;
v_y_hat         = sim_output.V_y_hat.signals.values;

time            = sim_output.eta.time;

%% Run plot script
plot_simulation_dp