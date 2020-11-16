clear all; close all; clc;

%% Constants for the matrices
% Azimuth Angles
alpha1  = 35; %deg
alpha2  = -35; %deg
alpha3  = 35; %deg
alpha4  = -35; %deg

% Convert angles to rad
alpha1  = deg2rad(alpha1);
alpha2  = deg2rad(alpha2);
alpha3  = deg2rad(alpha3);
alpha4  = deg2rad(alpha4);

%Surge position relative to the CO
l_x1    = 0.202; %metres
l_x2    = 0.202; %metres
l_x3    = -0.265; %metres
l_x4    = -0.265; %metres

% Sway position relative to the CO
l_y1    = -0.216; %metres
l_y2    = 0.216; %metres
l_y3    = 0.195; %metres
l_y4    = -0.195; %metres

% Rigid body mass
m_RB11  = 90; %kg
m_RB22  = 90; %kg
m_RB33  = 13; %kg

% Added mass
m_A11   = 54; %kg
m_A22   = 72; %kg
m_A33   = 5.2; %kg

% Damping coefficients
d11     = 250;
d22     = 200;
d33     = 15;

%% Definitions
b11     = cos(alpha1);
b12     = cos(alpha2);
b13     = cos(alpha3);
b14     = cos(alpha4);

b21     = sin(alpha1);
b22     = sin(alpha2);
b23     = sin(alpha3);
b24     = sin(alpha4);

b31     = l_x1*sin(alpha1)-l_y1*cos(alpha1);
b32     = l_x2*sin(alpha2)-l_y2*cos(alpha2);
b33     = l_x3*sin(alpha3)-l_y3*cos(alpha3);
b34     = l_x4*sin(alpha4)-l_y4*cos(alpha4);



%% System matrices for model equation
% Mass matrices
M_RB    = diag([m_RB11 m_RB22 m_RB33]);
M_RBinv = inv(M_RB);
M_A     = diag([m_A11 m_A22 m_A33]);
% Damping matrix
D       = diag([d11 d22 d33]);
% Actuator matrix or thrust allocation matrix
B       = [b11 b12 b13 b14;
            b21 b22 b23 b24;
            b31 b32 b33 b34];

% Coriolis-centripetal matrices
C_A     = @(ny_vec)([0 0 -m_A22*ny_vec(2);
                        0 0 m_A11*ny_vec(1);
                        m_A22*ny_vec(2) -m_A11*ny_vec(1) 0]);
                    
C_RB    = @(ny_vec)([0 0 -m_RB22*ny_vec(2);
                        0 0 m_RB11*ny_vec(1);
                        m_RB22*ny_vec(2) -m_RB11*ny_vec(1) 0]);

% Principle rotation matrix about z-axis
R_z     = @(psi)([cos(psi) -sin(psi) 0;
                  sin(psi) cos(psi) 0;
                  0 0 1]);
              
%% Initial values
% 3 DoF position in NED frame
x0      = 0; %m
y0      = 0; %m
psi0    = 0; %deg
psi0    = deg2rad(psi0); %rad

eta0    = [x0; y0; psi0];

% 3 DoF velocities in body frame
u0      = 0; % m/s
v0      = 0; % m/s
r0      = 0; % rad/s

ny0     = [u0; v0; r0];

% Adaptive initial values
theta_u_hat_init = [0; 0];
theta_v_hat_init = [0; 0];
theta_r_hat_init = [0; 0; 0; 0; 0];