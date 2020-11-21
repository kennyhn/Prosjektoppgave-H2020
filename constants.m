clear all; close all; clc;

%% Constants for the matrices
% ROV dimension
V_x         = 1;
V_y         = 0;

ROV_length  = 0.9;
ROV_width   = 0.65;
ROV_height  = 0.5;

% Water density
rho = 1000;

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

% These values are found through the master thesis of Herman Amundsen
% Rigid body mass
m_RB11  = 90; %kg
m_RB22  = 90; %kg
m_RB33  = 90; %kg
m_RB44  = 10; %kg m^2
m_RB55  = 15; % kg m^2
m_RB66  = 13; % kg m^2

% Added mass
m_A11   = 54; %kg
m_A22   = 72; %kg
m_A33   = 360; %kg
m_A44   = 11;  % kg
m_A55   = 43.5; %kg
m_A66   = 5.2; %kg

% Total mass = Added mass + Rigid body mass
m_11    = m_RB11+m_A11;
m_22    = m_RB22+m_A22;
m_33    = m_RB33+m_A33;
m_44    = m_RB44+m_A44;
m_55    = m_RB55+m_A55;
m_66    = m_RB66+m_A66;

% Damping coefficients
d_11    = 250;
d_22    = 200;
d_33    = 175;
d_44    = 20;
d_55    = 20;
d_66    = 15;

%% Definitions
b11     = cos(alpha1);
b12     = cos(alpha2);
b13     = cos(alpha3);
b14     = cos(alpha4);
b15     = 0;
b16     = 0;

b21     = sin(alpha1);
b22     = sin(alpha2);
b23     = sin(alpha3);
b24     = sin(alpha4);
b25     = 0;
b26     = 0;

b31     = 0;
b32     = 0;
b33     = 0;
b34     = 0;
b35     = 1;
b36     = 1;

b41     = 0;
b42     = 0;
b43     = 0;
b44     = 0;
b45     = 0;
b46     = 0;

b51     = 0;
b52     = 0;
b53     = 0;
b54     = 0;
b55     = 0.216;
b56     = -0.216;

b61     = l_x1*sin(alpha1)-l_y1*cos(alpha1);
b62     = l_x2*sin(alpha2)-l_y2*cos(alpha2);
b63     = l_x3*sin(alpha3)-l_y3*cos(alpha3);
b64     = l_x4*sin(alpha4)-l_y4*cos(alpha4);
b65     = 0;
b66     = 0;
%% System matrices for model equation
% Mass matrices
M_RB    = diag([m_RB11 m_RB22 m_RB33 m_RB44 m_RB55 m_RB66]);
M_A     = diag([m_A11 m_A22 m_A33 m_A44 m_A55 m_A66]);
M       = M_RB + M_A;
M_inv   = inv(M);

% Damping matrix
D       = diag([d_11 d_22 d_33 d_44 d_55 d_66]);
% Actuator matrix or thrust allocation matrix
B       = [b11 b12 b13 b14 b15 b16;
            b21 b22 b23 b24 b25 b26;
            b31 b32 b33 b34 b35 b36;
            b41 b42 b43 b44 b45 b46;
            b51 b52 b53 b54 b55 b56;
            b61 b62 b63 b64 b65 b66];


B_law       = [b11 b12 b13 b14 b15 b16;
            b21 b22 b23 b24 b25 b26;
            b31 b32 b33 b34 b35 b36;
            b61 b62 b63 b64 b65 b66];

B_psinv_law = B_law'/(B_law*B_law');



% Principle rotation matrix about z-axis
R_z     = @(psi)([cos(psi) -sin(psi) 0 0;
                  sin(psi) cos(psi) 0 0;
                  0 0 1 0;
                  0 0 0 1]);
              
%% Initial values
% 4 DoF position in inertial frame
x0      = 0; %m
y0      = 0; %m
z0      = 0; %m
phi0    = 0; %deg
theta0  = 0; %deg
psi0    = 0; %deg
phi0    = deg2rad(phi0);   %rad
theta0  = deg2rad(theta0); %rad
psi0    = deg2rad(psi0);   %rad

eta0    = [x0; y0; z0; phi0; theta0; psi0];

% 4 DoF velocities in body frame
u0      = 0; % m/s
v0      = 0; % m/s
w0      = 0; % m/s
p0      = 0; % rad/s
q0      = 0; % rad/s
r0      = 0; % rad/s

ny0     = [u0; v0; w0; p0; q0; r0];

% Thruster saturations
% Controller saturations
F_max       = 120; % N, maximum thruster force

torque_co   = (cross([l_x1; l_y1; 0], [F_max*cos(alpha1); F_max*sin(alpha1); 0]) +... %Nm
            cross([l_x2; l_y2; 0], [-F_max*cos(alpha2); -F_max*sin(alpha2); 0]) +...
            cross([l_x3; l_y3; 0], [-F_max*cos(alpha3); -F_max*sin(alpha3); 0]) +...
            cross([l_x4; l_y4; 0], [F_max*cos(alpha4); F_max*sin(alpha4); 0]));

tau_u_max   = (F_max*cos(alpha1)+F_max*cos(alpha2)+F_max*cos(alpha3)+F_max*cos(alpha4)); % N
tau_v_max   = (F_max*sin(alpha1)-F_max*sin(alpha2)+F_max*sin(alpha3)-F_max*sin(alpha4)); % N
tau_w_max   = 2*F_max; %N
tau_r_max   = sqrt(torque_co'*torque_co);