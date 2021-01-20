%% Controller gains
% Velocity controller gains
pole_u  = 6.736111111;
k_p_u   = m_11*pole_u-d_11;
k_i_u   = k_p_u/1000;
pole_v  = 6.234567901;
k_p_v   = m_22*pole_v-d_22;
k_i_v   = k_p_v/1000;

% Depth controller gains
zeta_d_heave = 0.5; % Critical damping
wb_d_heave   = 0.35; % desired bandwidth on heave
wn_heave     = wb_d_heave/sqrt(1-2*zeta_d_heave^2+sqrt(4*zeta_d_heave^4-4*zeta_d_heave^2+2));

k_p_z        = wn_heave^2*m_33;
k_i_z        = wn_heave/50*k_p_z;

% Heading controller gains
% Might need to do a new analysis on this tuning
zeta_d_psi   = 0.498236818; % Critical damping
wb_d_psi     = 2.039353174; % desired bandwidth on heading
wn_psi       = wb_d_psi/sqrt(1-2*zeta_d_psi^2+sqrt(4*zeta_d_psi^4-4*zeta_d_psi^2+2));

k_p_psi      = wn_psi^2*m_66;
k_d_psi      = 2*m_66*zeta_d_psi*wn_psi-d_66;
k_i_psi      = wn_psi/10*k_p_psi;

%% Adaptive gains
gamma_u     = 0.1;
gamma_v     = 0.1;
gamma_r     = 0.1;
lambda      = 1;