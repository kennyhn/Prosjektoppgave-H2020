% Options
use_saved_file = 0; % 1 for true 0 for false

if use_saved_file == 1
   clear all;
   step_response  = 1; % 0 for guidance, 1 for step
   filename   = 'simulation_output/dp_model/dp_model';
   
   if step_response == 1
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

%% Parse out simulation results
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
v_x_hat         = sim_output.V_x_hat.signals.values;
v_y_hat         = sim_output.V_y_hat.signals.values;

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
if step_response == 0
    plot([y_start y_los], [x_start x_los], 'g-x', 'LineWidth', 2);
end
plot(y, x, 'b', 'LineWidth', 2);
plot(x_nodes,y_nodes, 'r','LineWidth',3)
ylim([-2*50 2*50]);
xlim([-2*50 5*50]);
hold off;
grid on;
if step_response == 0
    legend_name = {'Path', 'ROV', 'Fish farm'};
else
    legend_name = {'ROV', 'Fish farm'};
end
legend(legend_name, 'FontSize', 14);
title('The position of the ROV in NED', 'FontSize', 16);
xlabel('y-EAST [m]');
ylabel('x-NORTH [m]');

figure();
gcf();
hold on;
plot(time, crosstrack_e, 'b', 'LineWidth', 2);
hold off;
grid on;
title('The crosstrack error', 'FontSize', 16);
legend('$e(t)$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('time [s]');
ylabel('error [m]', 'Interpreter', 'latex');

figure();
gcf();
hold on;
plot(time, z, 'b', time, z_d, 'r--', 'LineWidth', 2);
hold off;
grid on;
title('The depth the ROV goes', 'FontSize', 16);
legend('$z$', '$z_d$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('time [s]');
ylabel('z [m]', 'Interpreter', 'latex');

figure();
gcf();
hold on;
plot(time, rad2deg(psi), 'b', time, rad2deg(psi_d), 'r--', 'LineWidth', 2);
hold off;
grid on;
title('The angle of the ROV in NED', 'FontSize', 16);
legend('$\psi$', '$\psi_d$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('t [s]');
ylabel('$\psi$ [deg]', 'Interpreter', 'latex');

figure();
gcf();
hold on;
plot(time, u, 'b', time, u_d, 'r--', 'LineWidth', 2);
hold off;
grid on;
legend('$u$', '$u_d$' ,'Interpreter', 'latex', 'FontSize', 14);
title('The velocity in x-direction of the ROV in NED', 'FontSize', 16);
xlabel('t [s]');
ylabel('u [m/s]');

figure();
gcf();
hold on;
plot(time, v, 'b', time, v_d, 'r--', 'LineWidth', 2);
hold off;
grid on;
legend('$v$','$v_d$' ,'Interpreter', 'latex', 'FontSize', 14);
title('The velocity in y-direction of the ROV in NED', 'FontSize', 16);
xlabel('t [s]');
ylabel('v [m/s]');

figure();
gcf();
subplot(2, 2, 1);
hold on;
plot(time, tau_u_sat, 'b', 'LineWidth', 2);
plot(time, tau_u_unsat, 'r--', 'LineWidth', 2);
hold off;
grid on;
legend('$\tau_u$', '$\tau_{u_{unsat}}$', 'Interpreter' ,'latex');
title('Controller input surge');
xlabel('t [s]');
ylabel('$\tau_u$ [N]', 'Interpreter', 'latex');

subplot(2, 2, 2);
hold on;
plot(time, tau_v_sat, 'b', 'LineWidth', 2);
plot(time, tau_v_unsat, 'r--', 'LineWidth', 2);
hold off;
grid on;
legend('$\tau_v$', '$\tau_{v_{unsat}}$', 'Interpreter' ,'latex');
title('Controller input sway');
xlabel('t [s]');
ylabel('$\tau_v$ [N]', 'Interpreter', 'latex');

subplot(2, 2, 3);
hold on;
plot(time, tau_w_sat, 'b', 'LineWidth', 2);
plot(time, tau_w_unsat, 'r--', 'LineWidth', 2);
hold off;
grid on;
legend('$\tau_w$', '$\tau_{w_{unsat}}$', 'Interpreter', 'latex');
title('Controller input heave');
xlabel('t [s]');
ylabel('$\tau_w$ [Nm]', 'Interpreter', 'latex');

subplot(2, 2, 4);
hold on;
plot(time, tau_r_sat, 'b', 'LineWidth', 2);
plot(time, tau_r_unsat, 'r--', 'LineWidth', 2);
hold off;
grid on;
legend( '$\tau_r$', '$\tau_{r_{unsat}}$','Interpreter' ,'latex');
title('Controller input heading');
xlabel('t [s]');
ylabel('$\tau_r$ [Nm]', 'Interpreter', 'latex');
hold off;


figure();
gcf();
subplot(2, 1, 1);
hold on
plot(time, v_x_hat, 'b', 'LineWidth', 2);
plot(time, Vx, 'r--', 'LineWidth', 2);
hold off;
grid on;
legend('$\hat{V}_{x}$', '$V_x$', 'Interpreter' ,'latex');
title('Current estimate in x-direction, surge');
xlabel('t [s]');
ylabel('$V_x$ [m/s]', 'Interpreter', 'latex');

subplot(2, 1, 2);
hold on;
plot(time, v_y_hat, 'b', 'LineWidth', 2);
plot(time, Vy, 'r--', 'LineWidth', 2);
grid on;
legend('$\hat{V}_{y}$', '$V_y$', 'Interpreter' ,'latex');
title('Current estimate in y-direction, surge');
xlabel('t [s]');
ylabel('$V_y$ [m/s]', 'Interpreter', 'latex');