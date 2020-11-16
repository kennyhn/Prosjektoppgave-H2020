clc; close all;
%% Create all constants
constants

%% Disturbance
% Current (disturbance). Constant
Vx      = 0; % m/s in x-axis NED frame
Vy      = 0.2; % m/s in y-axis NED frame
Vc      = [Vx; Vy; 0];

ny_c0   = R_z(ny0(3))'*Vc;
% Controller saturations
F_max       = 120; % N, maximum thruster force

torque_co   = (cross([l_x1; l_y1; 0], [F_max*cos(alpha1); F_max*sin(alpha1); 0]) +... %Nm
            cross([l_x2; l_y2; 0], [-F_max*cos(alpha2); -F_max*sin(alpha2); 0]) +...
            cross([l_x3; l_y3; 0], [-F_max*cos(alpha3); -F_max*sin(alpha3); 0]) +...
            cross([l_x4; l_y4; 0], [F_max*cos(alpha4); F_max*sin(alpha4); 0]));

tau_u_max   = (F_max*cos(alpha1)+F_max*cos(alpha2)+F_max*cos(alpha3)+F_max*cos(alpha4))/(m_A11+m_RB11); % N
tau_v_max   = (F_max*sin(alpha1)-F_max*sin(alpha2)+F_max*sin(alpha3)-F_max*sin(alpha4))/(m_A22+m_RB22); % N
tau_r_max   = sqrt(torque_co'*torque_co)/(m_A33+m_RB33);

%% Controller gains
% PID gains

k_p_u = tau_u_max/100;
k_i_u = 0.1;
k_d_u = (2*1*k_p_u*1-d11/(m_A11+m_RB11));

k_p_v = tau_v_max/100;
k_i_v = 0.1;
k_d_v = (2*1*k_p_v*1-d22/(m_A22+m_RB22));

k_p_psi = tau_r_max/deg2rad(45);
k_d_psi = (2*1*k_p_psi*1-d33/(m_A33+m_RB33));



%Feedback linearizing gains
%{
% Surge controller
k_u         = 1; 

gamma_u     = 0.1;

% Sway controller
k_v         = 10;

gamma_v     = 0.1;

% Yaw controller
lambda      = 1;
k_psi       = 10;
k_r         = 10;

gamma_r     = 0.1;
%}


%% Simulation parameters
t_sim = 100; %s

%% Run simulation
sim_output = sim('simulering_ROV_feedback_linearizing.slx');

%% Parse out results
% For å sjekke at alt er likt uavhengig om man integrerer eller regner ut
% direkte
ny_c_direkt1    = sim_output.ny_c_direkt1;
ny_c_direkt2    = sim_output.ny_c_direkt2;
ny_c_direkt3    = sim_output.ny_c_direkt3;
ny_c_integrert1 = sim_output.ny_c_integrert1;
ny_c_integrert2 = sim_output.ny_c_integrert2;
ny_c_integrert3 = sim_output.ny_c_integrert3;

% Tilstander
x               = sim_output.x;
y               = sim_output.y;
psi             = sim_output.psi;
u               = sim_output.u;
v               = sim_output.v;
r               = sim_output.r;

% Referanser
u_d_out         = sim_output.u_d;
v_d_out         = sim_output.v_d;
psi_d_out       = sim_output.psi_d;


% Control inputs
tau_u_unsat     = sim_output.tau_u_unsat;
tau_v_unsat     = sim_output.tau_v_unsat;
tau_r_unsat     = sim_output.tau_r_unsat;

tau_u_sat       = sim_output.tau_u_sat;
tau_v_sat       = sim_output.tau_v_sat;
tau_r_sat       = sim_output.tau_r_sat;

% Parameterestimater for forstyrrelsene

%{
v_x_theta_u     = sim_output.theta_u_hat1;
v_y_theta_u     = sim_output.theta_u_hat2;

v_x_theta_v     = sim_output.theta_v_hat1;
v_y_theta_v     = sim_output.theta_v_hat2;
%}

time            = x.time;

%% Run plot script
plot_simulation