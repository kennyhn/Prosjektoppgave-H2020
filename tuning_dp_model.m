%% Controller options
int_effect_psi      = 1; % Turn on/off integral effect in heading controller
damping_estimation  = 1; % Turn on/off linear damping estimation on the controller or use linear values (assumed no nonlinear terms)

%% Controller gains
% Velocity controller gains
pole_u          = 6.736111111; % place pole at -value. Tune here
k_p_u           = pole_u-d_11/m_11;
pole_v          = 6.234567901; % place pole at -value. Tune here
k_p_v           = pole_v-d_22/m_22;

% Depth controller gains
zeta_d_heave    = 0.5; % Critical damping
wb_d_heave      = 0.35; % desired bandwidth on heave
wn_heave        = wb_d_heave/sqrt(1-2*zeta_d_heave^2+sqrt(4*zeta_d_heave^4-4*zeta_d_heave^2+2));

k_p_z           = wn_heave^2*m_33; 
k_d_z           = 2*m_33*zeta_d_heave*wn_heave-d_33;
k_i_z           = wn_heave/50*k_p_z;

% Heading controller gains
zeta_d_psi      = 0.498236818; % damping ratios
wb_d_psi        = 2.039353174; % desired bandwidth on heading
wn_psi          = wb_d_psi/sqrt(1-2*zeta_d_psi^2+sqrt(4*zeta_d_psi^4-4*zeta_d_psi^2+2));

k_p_psi         = wn_psi^2*m_66;
k_d_psi         = 2*m_66*zeta_d_psi*wn_psi-d_66;
k_i_psi         = wn_psi/10*k_p_psi;

%% Adaptive controller gains
gamma1          = 50;
gamma2          = 35;

%% Run simulation
sim_output = sim('simulering_ROV_DP_model.slx');
