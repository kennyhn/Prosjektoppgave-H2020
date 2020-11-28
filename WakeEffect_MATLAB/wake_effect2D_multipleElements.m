%**********************************************************************
%    Subroutine wake_effect2D_multipleElements
%    
%    Purpose: Top routine for calculaton of velocities through
%         and around a circular net.  
%    Method:  Mean turbulent wake velocites are calculated using 
%         linearized turbulent wake equations at individual cylinders 
%         that resemples the twines. Velocities in front of cage and  
%         outside of cage are calculated using lagally's theorem.  
%    Parameters
%    Input:
%    nu         - Kinematic viscosity
%    U_inf      - Uniform ambient velocity
%    Sn         - Solidity ratio
%    D_net      - Net diameter
%    d_tw       - Twine diameter
%    N_elements - Number of screen elements  
%    dx, dy     - Flow field spacial step length
%    xmin, xmax - Flow field end coordinates in x-direction 
%    ymin, ymax - Flow field end coordinates in y-direction  
%
%    Programmed by: Lars Haug
%    Date: September 2020
%*********************************************************************

clear all
close all
clc
addpath(genpath(pwd))
%% Run Settings
coefficientPlotting = true;
savefigure = true;
% Wake calculation method
% Sum of cylinders: 1
% 2D Plane wake: 2
method = 1;

% Paper comparisson
% TurnerEtAl: 1
% GanselEtAl: 2
PaperComparisson = 1;
%% System properties

nu = 1e-6;
rho_s = 1025;
if PaperComparisson == 1
    U_inf = 0.32;
    Sn = 0.3;
    D_net = 2.12;
    d_tw = 0.001;
    Re = U_inf*d_tw/(nu*(1-Sn));
    
elseif PaperComparisson == 2
    U_inf = 0.2;
    Sn = 0.25;
    D_net = 50;
    d_tw = 0.0015;
    Re = U_inf*d_tw/(nu*(1-Sn));
    
else
    Sn = 0.2;
    U_inf = 0.032;
    D_net = 50;
    d_tw = 0.0025;
    Re = U_inf*d_tw/(nu*(1-Sn));
end
% Find meshsize of net based on a square-shaped knotless mesh
p = [Sn -2*d_tw d_tw^2];
l = roots(p); 
[~,~,lambda] = find(l.*(l>d_tw)); % Choosing mesh size root larger than twine diameter
%Re_cyl = U_inf*D_net/nu;
%% Discretization of flow field
N_elements = 32;
dx = D_net/50; dy = D_net/50;
xmin = -2*D_net;
xmax = 5*D_net;
ymax = 2*D_net;

% Flow field grid
x_vec = xmin:dx:xmax;
y_vec = -ymax:dy:ymax;
%% Discretization of screen elements
[x_nodes, y_nodes, lj, theta, phi] = elementDiscretization(N_elements, D_net);
Ni = floor(lj/lambda);
N_cyl = Ni*N_elements;
[xi, yi, exi, eyi] = sourcePointDiscretization(x_nodes, y_nodes, N_elements, Ni);
%% Memory Allocation
C_N = zeros(1, N_elements);
C_T = zeros(1, N_elements);
C_D = zeros(1, N_elements);
C_L = zeros(1, N_elements);

%% Calculation of inflow velocity on downstream net section
for e = N_elements/2+1:N_elements
    [C_N(e), C_T(e), C_D(e), C_L(e), cd, cl] = reynoldsDependentDragCoefficients(Re, Sn, phi(e));
end
U_dnet = U_inf*(1-0.46*cd);
if U_dnet < 0
    disp('Warning: The inflow velocity to the net is less than zero')
    pause
end
Re_dnet = U_dnet.*d_tw./(nu*(1-Sn));
%% Coefficients
cd_loland = 0.04 + (-0.04+0.33*Sn+6.54*Sn^2-4.88*Sn^3);
for e = 1:N_elements/2
    [C_N(e), C_T(e), C_D(e), C_L(e)] = reynoldsDependentDragCoefficients(Re_dnet, Sn, phi(e));
    C_Di(1,(e-1)*Ni+1:e*Ni) = C_D(e)/Sn;
end
%% Flow Around Cage and wake (Lagally's theorem and 2D plane wake)
[u,v, PHI, PSI, u1_pWake] = wakeFlow(U_dnet, D_net, N_elements, U_inf, rho_s, C_D, lj, x_vec, y_vec, xi, yi, C_Di, d_tw, Ni, cd, dy, ymax, method);
%% Post Processing
disp(['Plotting...'])
% Coefficent plotting
if coefficientPlotting == true
    coefficientPlot(N_elements,C_N,C_D,C_T,C_L,phi,cl,cd)
end
topNodeIndx = int64(1/2*N_elements+1);
boundary(:,1) = 0:dx:xmax;
A = transpose( 0.95 * sqrt(C_D(topNodeIndx-1) / Sn * d_tw / (boundary(:,1) - x_nodes(topNodeIndx)) ) );
A = A+(A<=0.01)*0.01;
B = 0.0888 * C_D(topNodeIndx-1) * d_tw * (boundary(:,1) - x_nodes(topNodeIndx));
boundary(:,2) = y_nodes(topNodeIndx) + sqrt(B .* (log(A) - log(0.01)));

% Flow Field Plotting
flowPlot(xmin, xmax, ymax, D_net, dx, dy, x_vec, y_vec, x_nodes, y_nodes, U_dnet, U_inf, u, v, PSI, cd, cd_loland, PaperComparisson,boundary,savefigure)
disp(['Finished'])