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
    t_sim = 1450; %s
end

%% Run simulation
sim_output = sim('simulering_ROV_FL_controller.slx');