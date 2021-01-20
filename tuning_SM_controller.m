%% Controller gains
% Velocity controller gains
lambda_u    = 0.5; %10;
phi_u       = 1; %5;
k_d_u       = 5;
k_s_u       = 500; %3;

lambda_v    = 0.5; %10;
phi_v       = 1; %5;
k_d_v       = 3; %5;
k_s_v       = 400; %3;

% Depth controller gains
lambda_z    = 3; %10;
phi_z       = 0.5; %1;
k_d_z       = 30; %1; 
k_s_z       = 10; %1;

% Heading controller gains
lambda_psi  = 1; %35;
phi_psi     = 0.5; %deg2rad(10);
k_d_psi     = 10; %20;
k_s_psi     = 50; %25;

%% Adaptive gains
gamma_u     = 0.1;
gamma_v     = 0.1;
gamma_r     = 0.1;
lambda      = 1;