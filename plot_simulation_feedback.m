%% Plot the result
figure();
gcf();
hold on;
plot(y, x, 'b');
hold off;
grid on;
legend('ROV');
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
plot(time, rad2deg(psi), 'b', time, rad2deg(psi_d), 'b--');
hold off;
grid on;
title('The angle of the ROV in NED');
legend('$\psi$', 'interpreter', 'latex');
xlabel('t [s]');
ylabel('$\psi$ [deg]', 'interpreter', 'latex');

figure();
gcf();
hold on;
plot(time, u, 'b');
hold off;
grid on;
legend('u', 'Interpreter', 'latex');
title('The velocity in x-direction of the ROV in NED');
xlabel('t [s]');
ylabel('u [m/s]');

figure();
gcf();
hold on;
plot(time, v , 'b');
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
plot(time, tau_u_sat, 'b--');
plot(time, tau_u_unsat, 'b');
hold off;
grid on;
legend('$\tau_u$', '$\tau_{u_{unsat}}$', 'Interpreter' ,'latex');
title('Controller input surge');
xlabel('t [s]');
ylabel('$\tau_u$ [N]', 'Interpreter', 'latex');

subplot(2, 2, 2);
hold on;
plot(time, tau_v_sat, 'b--');
plot(time, tau_v_unsat, 'b');
hold off;
grid on;
legend('$\tau_v$', '$\tau_{v_{unsat}}$', 'Interpreter' ,'latex');
title('Controller input sway');
xlabel('t [s]');
ylabel('$\tau_v$ [N]', 'Interpreter', 'latex');

subplot(2, 2, 3);
hold on;
plot(time, tau_w_sat, 'b--');
plot(time, tau_w_unsat, 'b');
hold off;
grid on;
legend('$\tau_w$', '$\tau_{w_{unsat}}$', 'Interpreter', 'latex');
title('Controller input heave');
xlabel('t [s]');
ylabel('$\tau_w$ [Nm]', 'Interpreter', 'latex');

subplot(2, 2, 4);
hold on;
plot(time, tau_r_sat, 'b--');
plot(time, tau_r_unsat, 'b');
hold off;
grid on;
legend('$\tau_r$', '$\tau_{r_{unsat}}$', 'Interpreter' ,'latex');
title('Controller input heading');
xlabel('t [s]');
ylabel('$\tau_r$ [Nm]', 'Interpreter', 'latex');
hold off;


figure();
gcf();
subplot(2, 2, 1);
hold on
plot(time, Vx*ones(size(time)), 'b--');
plot(time, v_xu_hat, 'b');
hold off;
grid on;
legend('$V_x$', '$\hat{V}_{xu}$', 'Interpreter' ,'latex');
title('Current estimate in x-direction, surge');
xlabel('t [s]');
ylabel('$V_x$ [m/s]', 'Interpreter', 'latex');

subplot(2, 2, 2);
hold on;
plot(time, Vy*ones(size(time)), 'b--');
plot(time, v_yu_hat, 'b');
grid on;
legend('$V_y$', '$\hat{V}_{yu}$', 'Interpreter' ,'latex');
title('Current estimate in y-direction, surge');
xlabel('t [s]');
ylabel('$V_y$ [m/s]', 'Interpreter', 'latex');

subplot(2, 2, 3);
hold on;
plot(time, Vx*ones(size(time)), 'b--');
plot(time, v_xv_hat, 'b');
grid on;
legend('$V_x$', '$\hat{V}_{xv}$', 'Interpreter' ,'latex');
title('Current estimate in x-direction, sway');
xlabel('t [s]');
ylabel('$V_x$ [m/s]', 'Interpreter', 'latex');

subplot(2, 2, 4);
hold on;
plot(time, Vy*ones(size(time)), 'b--');
plot(time, v_yv_hat, 'b');
grid on;
legend('$V_y$', '$\hat{V}_{yv}$', 'Interpreter' ,'latex');
title('Current estimate in y-direction, sway');
xlabel('t [s]');
ylabel('$V_y$ [m/s]', 'Interpreter', 'latex');
hold off;


figure();
gcf();
subplot(2, 2, 1);
hold on
plot(time, Vx*ones(size(time)), 'b--');
plot(time, v_xr_hat, 'b');
hold off;
grid on;
legend('$V_x$', '$\hat{V}_{xr}$', 'Interpreter' ,'latex');
title('Current estimate in x-direction, yaw');
xlabel('t [s]');
ylabel('$V_x$ [m/s]', 'Interpreter', 'latex');

subplot(2, 2, 2);
hold on;
plot(time, Vy*ones(size(time)), 'b--');
plot(time, v_yr_hat, 'b');
grid on;
legend('$V_y$', '$\hat{V}_{yr}$', 'Interpreter' ,'latex');
title('Current estimate in y-direction, yaw');
xlabel('t [s]');
ylabel('$V_y$ [m/s]', 'Interpreter', 'latex');

subplot(2, 2, 3);
hold on;
plot(time, Vx*Vx*ones(size(time)), 'b--');
plot(time, v_xr_hat_sq, 'b');
grid on;
legend('$V_x^2$', '$\hat{V}_{xr}^2$', 'Interpreter' ,'latex');
title('Current estimate in x-direction, yaw');
xlabel('t [s]');
ylabel('$V_x^2$ [m/s]', 'Interpreter', 'latex');

subplot(2, 2, 4);
hold on;
plot(time, Vy*Vy*ones(size(time)), 'b--');
plot(time, v_yr_hat_sq, 'b');
grid on;
legend('$V_y^2$', '$\hat{V}_{yr}^2$', 'Interpreter' ,'latex');
title('Current estimate in y-direction, yaw');
xlabel('t [s]');
ylabel('$V_y^2$ [m/s]', 'Interpreter', 'latex');
hold off;

figure()
gcf();
hold on;
plot(time, Vx*Vy*ones(size(time)), 'b--');
plot(time, v_xyr_hat, 'b');
legend('$V_{xy}$', '$\hat{V}_{xyr}$', 'Interpreter' ,'latex');
grid on;
title('Current estimate xy, yaw');
xlabel('t [s]');
ylabel('$V_{xy}$ [m/s]', 'Interpreter', 'latex');
hold off;