clc; close all;
%% Create all constants
constants_model_coriolis

%% Disturbance
% Current (disturbance). Constant
V_x     = 0.25; % m/s in x-axis inertial frame
V_y     = -0.6; % m/s in y-axis inertial frame
g_z     = 9.81;


%% Controller gains
% Velocity controller gains
pole_u      = 5; % place pole at -value. Tune here
k_p_u       = m_11*pole_u-d_11;
gamma_u     = 0.1;
pole_v      = 5; % place pole at -value. Tune here
k_p_v       = m_22*pole_v-d_22;
gamma_v     = 0.1;


% Depth controller gains
zeta_d_heave = 1; % Critical damping
wb_d_heave   = 1; % desired bandwidth on heave
wn_heave     = wb_d_heave/sqrt(1-2*zeta_d_heave^2+sqrt(4*zeta_d_heave^4-4*zeta_d_heave^2+2));

k_p_w        = wn_heave^2*m_33; 
k_d_w        = 2*m_33*zeta_d_heave*wn_heave-d_33;
k_i_w        = wn_heave/10*k_p_w;

% Heading controller gains
lambda       = d_44; % tuning system damping must be greater than 0
zeta_d_psi   = 1; % Critical damping
wb_d_psi     = 0.8; % desired bandwidth on heading
wn_psi       = wb_d_psi/sqrt(1-2*zeta_d_psi^2+sqrt(4*zeta_d_psi^4-4*zeta_d_psi^2+2));

k_p_r        = 2*zeta_d_psi*wn_psi-lambda;
k_p_psi      = wn_psi^2-lambda*k_p_r;
gamma_r      = 0.1;


%% Simulation parameters
t_sim = 150; %s

%% Run simulation
sim_output = sim('simulering_ROV_FL_controller.slx');

%% Parse out results
nu              = sim_output.nu.signals.values;
eta             = sim_output.eta.signals.values;

tau_unsat       = sim_output.tau_unsat.signals.values;
tau_sat         = sim_output.tau_sat.signals.values;

psi_d           = sim_output.psi_d.signals.values;

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