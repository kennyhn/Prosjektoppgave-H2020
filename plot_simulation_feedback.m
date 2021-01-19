% Options
use_saved_file = 0; % 1 for true 0 for false


if use_saved_file == 1
   clear; close all
   % options
   controller     = 0; % 0 for FL, 1 for PID, 2 for SM
   step_response  = 0; % 0 for guidance, 1 for step
   filename   = 'simulation_output/';
   
   if controller == 0
       filename = strcat(filename, 'FL_controller/FL_controller');
   elseif controller == 1
       filename = strcat(filename, 'PID_controller/PID_controller');
   elseif controller == 2
       filename = strcat(filename, 'SM_controller/SM_controller');
   end
   
   if step_response == 0
       % Guidance law parameters
       Delta   = 25; % Lookahead distance
       x_start = -2*50;
       y_start = 1*80;
       x_los   = 2*50;
       y_los   = 1*80;

       filename = strcat(filename, '_guidance.mat');
   else  
       filename = strcat(filename, '_step.mat');
   end
   
   sim_output = load(filename).sim_output;
end

%% Parse out results from simulation
nu              = sim_output.nu.signals.values;
eta             = sim_output.eta.signals.values;

tau_unsat       = sim_output.tau_unsat.signals.values;
tau_sat         = sim_output.tau_sat.signals.values;

psi_d           = sim_output.psi_d.signals.values;
u_d             = sim_output.u_d.signals.values;
v_d             = sim_output.v_d.signals.values;
z_d             = sim_output.z_d.signals.values;

Vc              = sim_output.disturbance.signals.values;
Vx              = Vc(:, 1);
Vy              = Vc(:, 2);

crosstrack_e    = sim_output.crosstrack_error.signals.values;

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

%% Create the circle representing the fish farm
N_elements  = 1000; % Size of points
N_nodes     = N_elements +1;
D_net       = 50; % Diameter of the fish farm net
dtheta      = 2*pi/N_elements;
theta       = (0:dtheta:N_elements*dtheta)-N_elements/2*dtheta/2;
[x_nodes, y_nodes] = pol2cart(theta, ones(1,N_nodes)*D_net/2);

%% Plot the result
figure();
gcf();
hold on;
plot(y, x, 'b', 'LineWidth', 2);
if step_response == 0
    plot(waypoints(:, 2), waypoints(:, 1), 'g--x', 'LineWidth', 2);
end
plot(x_nodes,y_nodes, 'r','LineWidth',3)
hold off;
ylim([-2*50 2*50]);
xlim([-2*50 5*50]);
grid on;
if step_response == 0
    legend_name = {'ROV','Path','Fish farm'};
else
    legend_name = {'ROV', 'Fish farm'};
end
legend(legend_name, 'FontSize', 14);
title('\textbf{The position of the ROV in NED}','Interpreter', 'latex', 'FontSize', 16);
xlabel('y-EAST [m]');
ylabel('x-NORTH [m]');
set(gca, 'FontSize', 14);

figure();
gcf();
hold on;
plot(time, crosstrack_e, 'b', 'LineWidth', 2);
hold off;
xlim([0 time(end)]);
grid on;
title('\textbf{The cross-track error}','Interpreter', 'latex', 'FontSize', 16);
legend('$e(t)$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('time [s]');
ylabel('error [m]', 'Interpreter', 'latex');
set(gca, 'FontSize', 14);

figure();
gcf();
hold on;
plot(time, z, 'b', time, z_d, 'r--', 'LineWidth', 2);
hold off;
xlim([0 time(end)]);
grid on;
title('\textbf{The depth of the ROV relative to NED}', 'Interpreter', 'latex', 'FontSize', 16);
legend('$z$','$z_d$', 'interpreter', 'latex', 'FontSize', 14);
xlabel('time [s]');
ylabel('z [m]', 'interpreter', 'latex');
set(gca, 'FontSize', 14);

figure();
gcf();
hold on;
plot(time, rad2deg(psi), 'b', time, rad2deg(psi_d), 'r--', 'LineWidth', 2);
hold off;
xlim([0 time(end)]);
grid on;
title('\textbf{The yaw angle of the ROV relative to NED}', 'Interpreter', 'latex', 'FontSize', 16);
legend('$\psi$', '$\psi_d$', 'interpreter', 'latex', 'FontSize', 14);
xlabel('t [s]');
ylabel('$\psi$ [deg]', 'interpreter', 'latex');
set(gca, 'FontSize', 14);

figure();
gcf();
subplot(2, 1, 1);
hold on;
plot(time, u, 'b', time, u_d, 'r--', 'LineWidth', 2);
hold off;
xlim([0 time(end)]);
grid on;
legend('$u$', '$u_d$', 'Interpreter', 'latex', 'FontSize', 14);
title('\textbf{The surge velocity of the ROV in body}', 'interpreter', 'latex', 'FontSize', 16);
xlabel('t [s]');
ylabel('u [m/s]');
set(gca, 'FontSize', 14);
subplot(2, 1, 2);
hold on;
plot(time, v , 'b', time, v_d, 'r--', 'LineWidth', 2);
hold off;
xlim([0 time(end)]);
grid on;
legend('$v$', '$v_d$', 'Interpreter', 'latex', 'FontSize', 14);
title('\textbf{The sway velocity of the ROV in body}', 'interpreter', 'latex', 'FontSize', 16);
xlabel('t [s]');
ylabel('v [m/s]');
set(gca, 'FontSize', 14);

figure();
gcf();
subplot(2, 2, 1);
hold on;
plot(time, tau_u_unsat, 'r--', 'LineWidth', 2);
plot(time, tau_u_sat, 'b', 'LineWidth', 2);
hold off;
xlim([0 time(end)]);
grid on;
legend('$\tau_{u_{unsat}}$', '$\tau_u$', 'Interpreter', 'latex', 'FontSize', 14);
title('\textbf{Control input surge}', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('t [s]');
ylabel('$\tau_u$ [N]', 'Interpreter', 'latex');
set(gca, 'FontSize', 14);

subplot(2, 2, 2);
hold on;
plot(time, tau_v_unsat, 'r--', 'LineWidth', 2);
plot(time, tau_v_sat, 'b', 'LineWidth', 2);
hold off;
xlim([0 time(end)]);
grid on;
legend('$\tau_{v_{unsat}}$', '$\tau_v$', 'Interpreter' ,'latex', 'FontSize', 14);
title('\textbf{Control input sway}', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('t [s]');
ylabel('$\tau_v$ [N]', 'Interpreter', 'latex');
set(gca, 'FontSize', 14);

subplot(2, 2, 3);
hold on;
plot(time, tau_w_unsat, 'r--', 'LineWidth', 2);
plot(time, tau_w_sat, 'b', 'LineWidth', 2);
hold off;
xlim([0 time(end)]);
grid on;
legend('$\tau_{w_{unsat}}$', '$\tau_{w}$', 'Interpreter', 'latex', 'FontSize', 14);
title('\textbf{Control input heave}', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('t [s]');
ylabel('$\tau_w$ [Nm]', 'Interpreter', 'latex');
set(gca, 'FontSize', 14);

subplot(2, 2, 4);
hold on;
plot(time, tau_r_unsat, 'r--', 'LineWidth', 2);
plot(time, tau_r_sat, 'b', 'LineWidth', 2);
hold off;
xlim([0 time(end)]);
grid on;
legend('$\tau_{r_{unsat}}$', '$\tau_r$', 'Interpreter', 'latex', 'FontSize', 14);
title('\textbf{Control input heading}', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('t [s]');
ylabel('$\tau_r$ [Nm]', 'Interpreter', 'latex');
set(gca, 'FontSize', 14);
hold off;


figure();
gcf();
subplot(2, 2, 1);
hold on
plot(time, v_xu_hat, 'b', 'LineWidth', 2);
plot(time, Vx, 'r--', 'LineWidth', 2);
hold off;
xlim([0 time(end)]);
grid on;
legend('$\hat{V}_{xu}$', '$V_x$', 'Interpreter' ,'latex', 'FontSize', 14);
title('$\mathbf{u}$ \textbf{estimate in NED,} $\mathbf{\theta_u}$', 'Interpreter', 'Latex', 'FontSize', 16);
xlabel('t [s]');
ylabel('$V_x$ [m/s]', 'Interpreter', 'latex');
set(gca, 'FontSize', 14);

subplot(2, 2, 2);
hold on;
plot(time, v_yu_hat, 'b', 'LineWidth', 2);
plot(time, Vy, 'r--', 'LineWidth', 2);
hold off;
xlim([0 time(end)]);
grid on;
legend('$\hat{V}_{yu}$', '$V_y$', 'Interpreter' ,'latex', 'FontSize', 14);
title('$\mathbf{v}$ \textbf{estimate in NED,} $\mathbf{\theta_u}$', 'Interpreter', 'Latex', 'FontSize', 16);
xlabel('t [s]');
ylabel('$V_y$ [m/s]', 'Interpreter', 'latex');
set(gca, 'FontSize', 14);

subplot(2, 2, 3);
hold on;
plot(time, v_xv_hat, 'b', 'LineWidth', 2);
plot(time, Vx, 'r--', 'LineWidth', 2);
hold off;
xlim([0 time(end)]);
grid on;
legend('$\hat{V}_{xv}$', '$V_x$', 'Interpreter' ,'latex', 'FontSize', 14);
title('$\mathbf{u}$ \textbf{estimate in NED,} $\mathbf{\theta_v}$', 'Interpreter', 'Latex', 'FontSize', 16);
xlabel('t [s]');
ylabel('$V_x$ [m/s]', 'Interpreter', 'latex');
set(gca, 'FontSize', 14);

subplot(2, 2, 4);
hold on;
plot(time, v_yv_hat, 'b', 'LineWidth', 2);
plot(time, Vy, 'r--', 'LineWidth', 2);
hold off;
xlim([0 time(end)]);
grid on;
legend('$\hat{V}_{yv}$', '$V_y$', 'Interpreter' ,'latex', 'FontSize', 14);
title('$\mathbf{v}$ \textbf{estimate in NED,} $\mathbf{\theta_v}$', 'Interpreter', 'Latex', 'FontSize', 16);
xlabel('t [s]');
ylabel('$V_y$ [m/s]', 'Interpreter', 'latex');
set(gca, 'FontSize', 14);

figure();
gcf();
subplot(2, 2, 1);
hold on
plot(time, v_xr_hat, 'b', 'LineWidth', 2);
plot(time, Vx, 'r--', 'LineWidth', 2);
hold off;
xlim([0 time(end)]);
grid on;
legend('$\hat{V}_{xr}$', '$V_x$', 'Interpreter' ,'latex', 'FontSize', 14);
title('$\mathbf{u}$ \textbf{estimate in NED,} $\mathbf{\theta_r}$', 'Interpreter', 'Latex', 'FontSize', 16);
xlabel('t [s]');
ylabel('$V_x$ [m/s]', 'Interpreter', 'latex');
set(gca, 'FontSize', 14);

subplot(2, 2, 2);
hold on;
plot(time, v_yr_hat, 'b', 'LineWidth', 2);
plot(time, Vy, 'r--', 'LineWidth', 2);
hold off;
xlim([0 time(end)]);
grid on;
legend('$\hat{V}_{yr}$', '$V_y$', 'Interpreter' ,'latex', 'FontSize', 14);
title('$\mathbf{v}$ \textbf{estimate in NED,} $\mathbf{\theta_r}$', 'Interpreter', 'Latex', 'FontSize', 16);
xlabel('t [s]');
ylabel('$V_y$ [m/s]', 'Interpreter', 'latex');
set(gca, 'FontSize', 14);

subplot(2, 2, 3);
hold on;
plot(time, v_xr_hat_sq, 'b', 'LineWidth', 2);
plot(time, Vx.*Vx, 'r--', 'LineWidth', 2);
hold off;
xlim([0 time(end)]);
grid on;
legend('$\hat{V}_{xr}^2$', '$V_x^2$', 'Interpreter' ,'latex', 'FontSize', 14);
title('$\mathbf{u^2}$ \textbf{estimate in NED,} $\mathbf{\theta_r}$', 'Interpreter', 'Latex', 'FontSize', 16);
xlabel('t [s]');
ylabel('$V_x^2$ [m/s]', 'Interpreter', 'latex');
set(gca, 'FontSize', 14);

subplot(2, 2, 4);
hold on;
plot(time, v_yr_hat_sq, 'b', 'LineWidth', 2);
plot(time, Vy.*Vy, 'r--', 'LineWidth', 2);
hold off;
xlim([0 time(end)]);
grid on;
legend('$\hat{V}_{yr}^2$', '$V_y^2$', 'Interpreter' ,'Latex', 'FontSize', 14);
title('$\mathbf{v^2}$ \textbf{estimate in NED,} $\mathbf{\theta_r}$', 'Interpreter', 'Latex', 'FontSize', 16);
xlabel('t [s]');
ylabel('$V_y^2$ [m/s]', 'Interpreter', 'latex');
set(gca, 'FontSize', 14);

figure()
gcf();
hold on;
plot(time, v_xyr_hat, 'b', 'LineWidth', 2);
plot(time, Vx.*Vy, 'r--', 'LineWidth', 2);
hold off;
xlim([0 time(end)]);
legend('$\hat{V}_{xyr}$', '$V_{xy}$', 'Interpreter' ,'Latex', 'FontSize', 14);
grid on;
title('$\mathbf{uv}$ \textbf{estimate in NED,} $\mathbf{\theta_r}$', 'Interpreter', 'Latex', 'FontSize', 16);
xlabel('t [s]');
ylabel('$V_{xy}$ [m/s]', 'Interpreter', 'Latex');
set(gca, 'FontSize', 14);

[N, m] = size(u);


u_e = sumabs(u-u_d)/N
v_e = sumabs(v-v_d)/N
c_e = sumabs(crosstrack_e)/N