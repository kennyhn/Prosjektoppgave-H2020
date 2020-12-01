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
plot([y_start y_los], [x_start x_los], 'g-x');
plot(y, x, 'b');
plot(x_nodes,y_nodes, 'r','LineWidth',3)
ylim([-2*50 2*50]);
xlim([-2*50 5*50]);
hold off;
grid on;
legend('Path','ROV', 'Fish farm');
title('The position of the ROV in NED');
xlabel('y [m]');
ylabel('x [m]');

figure();
gcf();
hold on;
plot(time, z, 'b');
hold off;
grid on;
title('The depth the ROV goes');
legend('$z$', 'interpreter', 'latex');
xlabel('time [s]');
ylabel('z [m]', 'interpreter', 'latex');

figure();
gcf();
hold on;
plot(time, rad2deg(psi), 'b', time, rad2deg(psi_d), 'r--');
hold off;
grid on;
title('The angle of the ROV in NED');
legend('$\psi$', 'interpreter', 'latex');
xlabel('t [s]');
ylabel('$\psi$ [deg]', 'interpreter', 'latex');

figure();
gcf();
hold on;
plot(time, u, 'b', time, u_d, 'r--');
hold off;
grid on;
legend('u', 'Interpreter', 'latex');
title('The velocity in x-direction of the ROV in NED');
xlabel('t [s]');
ylabel('u [m/s]');

figure();
gcf();
hold on;
plot(time, v, 'b', time, v_d, 'r--');
hold off;
grid on;
legend('$v$', 'Interpreter', 'latex');
title('The velocity in y-direction of the ROV in NED');
xlabel('t [s]');
ylabel('v [m/s]');

figure();
gcf();
subplot(2, 2, 1);
hold on;
plot(time, tau_u_unsat, 'r--');
plot(time, tau_u_sat, 'b');
hold off;
grid on;
legend('$\tau_{u_{unsat}}$', '$\tau_u$', 'Interpreter' ,'latex');
title('Controller input surge');
xlabel('t [s]');
ylabel('$\tau_u$ [N]', 'Interpreter', 'latex');

subplot(2, 2, 2);
hold on;
plot(time, tau_v_unsat, 'r--');
plot(time, tau_v_sat, 'b');
hold off;
grid on;
legend('$\tau_{v_{unsat}}$', '$\tau_v$', 'Interpreter' ,'latex');
title('Controller input sway');
xlabel('t [s]');
ylabel('$\tau_v$ [N]', 'Interpreter', 'latex');

subplot(2, 2, 3);
hold on;
plot(time, tau_w_unsat, 'r--');
plot(time, tau_w_sat, 'b');
hold off;
grid on;
legend('$\tau_{w_{unsat}}$', '$\tau_w$', 'Interpreter', 'latex');
title('Controller input heave');
xlabel('t [s]');
ylabel('$\tau_w$ [Nm]', 'Interpreter', 'latex');

subplot(2, 2, 4);
hold on;
plot(time, tau_r_unsat, 'r--');
plot(time, tau_r_sat, 'b');
hold off;
grid on;
legend('$\tau_{r_{unsat}}$', '$\tau_r$', 'Interpreter' ,'latex');
title('Controller input heading');
xlabel('t [s]');
ylabel('$\tau_r$ [Nm]', 'Interpreter', 'latex');
hold off;


figure();
gcf();
subplot(2, 1, 1);
hold on
plot(time, Vx, 'r--');
plot(time, v_x_hat, 'b');
hold off;
grid on;
legend('$V_x$', '$\hat{V}_{x}$', 'Interpreter' ,'latex');
title('Current estimate in x-direction, surge');
xlabel('t [s]');
ylabel('$V_x$ [m/s]', 'Interpreter', 'latex');

subplot(2, 1, 2);
hold on;
plot(time, Vy, 'r--');
plot(time, v_y_hat, 'b');
grid on;
legend('$V_y$', '$\hat{V}_{y}$', 'Interpreter' ,'latex');
title('Current estimate in y-direction, surge');
xlabel('t [s]');
ylabel('$V_y$ [m/s]', 'Interpreter', 'latex');