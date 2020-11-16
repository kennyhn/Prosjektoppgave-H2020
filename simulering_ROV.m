clc; close all; 

%% Initialize the ROVs constants
constants

%% Simulation parameters
h = 0.01; % s
N = 1500; % steps
t = 0:h:(N-1)*h;

%% Time-varying variables
% Initialize the states
ny          = zeros(3, N);
ny_c        = zeros(3, N);
ny_r        = zeros(3, N);
eta         = zeros(3, N);

ny_dot      = zeros(3, 1);
ny_c_dot    = zeros(3, 1);
ny_r_dot    = zeros(3, 1);
eta_dot     = zeros(3, 1);

% Insert initial values
ny(:,1)     = ny0(:, 1);
ny_c(:, 1)  = R_z(ny(3, 1))'*Vc; % Convert the ocean current to body frame
ny_r(:, 1)  = ny(:, 1) - ny_c(:, 1);
eta(:, 1)   = eta0(:, 1);

% TBD: Make a more sophisticated controller input for f using control
% theory

% Controller inputs
f       = zeros(size(B, 2), N);

for i = 1:N-1
   % Calculate the change in states: eta_dot and ny_dot
   eta_dot = R_z(eta(3, i))*ny(:, i);
   
   ny_dot  = M_RBinv*(B*f(:,i)-C_RB(ny(:, i))*ny(:, i)-...
       M_A*ny_r_dot-C_A(ny_r(:, i))*ny_r(:, i)-D*ny_r(:,i));
   
   ny_c_dot= [ny(3, i)*ny_c(2,i); -ny(3, i)*ny_c(1, i); 0];
   
   % Update the states using eulers method
   eta(:, i+1)  = eta(:, i) + h*eta_dot;
   ny(:, i+1)   = ny(:, i) + h*ny_dot;
   
   ny_c(:, i+1) = R_z(eta(3, i+1))'*Vc;
   ny_r(:, i+1) = ny(:, i+1) - ny_c(:, i+1);
   
   ny_r_dot     = ny_dot-ny_c_dot;
   
   %% TBD: Controller
end


%% Plot the states
figure(1);
hold on;
plot(eta(1, :), eta(2, :), 'b');
hold off;
grid on;
legend('ROV');
title('The position of the ROV in NED');
xlabel('x [m]');
ylabel('y [m]');

figure(2);
hold on;
plot(t, rad2deg(eta(3, :)), 'b');
hold off;
grid on;
title('The angle of the ROV in NED');
legend('$\psi$', 'interpreter', 'latex');
xlabel('t [s]');
ylabel('$\psi$ [deg]', 'interpreter', 'latex');

figure(3);
hold on;
plot(t, ny(1, :), 'b');
hold off;
grid on;
legend('u');
title('The velocity in x-direction of the ROV in NED');
xlabel('t [s]');
ylabel('u [m/s]');

figure(4);
hold on;
plot(t, ny(2, :), 'b');
hold off;
grid on;
legend('v');
title('The velocity in y-direction of the ROV in NED');
xlabel('t [s]');
ylabel('v [m/s]');