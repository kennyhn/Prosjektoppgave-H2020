function [u, v, PHI, PSI] = wakeFlow(U_dStream, D_net, Ni, N_elements, U_inf, rho_s, C_D, lj, x_vec, y_vec, x_nodes, y_nodes)

error_param = 0.5;

% Memory allocation
PHI_0 = zeros(length(x_vec),length(y_vec));
PSI_0 = zeros(length(x_vec),length(y_vec));

PHI_source = zeros(length(x_vec),length(y_vec));
PSI_source = zeros(length(x_vec),length(y_vec));

u_source = zeros(length(x_vec),length(y_vec));
v_source = zeros(length(x_vec),length(y_vec));

u1_pWake = zeros(length(x_vec),length(y_vec));

U_elem = [U_dStream(1:Ni:end) ones(1,N_elements/2)*U_inf];
F_D = 0.5*rho_s.*U_elem.^2.*C_D*lj;

C_D_tot = 2*sum(F_D)/(rho_s*U_inf^2*D_net);
epsilon_tau = error_param*0.0222*C_D_tot*D_net*U_inf;

q_source = sum(F_D)/(rho_s*U_inf);

for k = 1:length(x_vec)
    x = x_vec(k);
    for i = 1:length(y_vec)
        y = y_vec(i);
        [th1, r1] = cart2pol(x,y);
        %Skip Calculations inside cage
        if x^2+y^2 < (D_net/2)^2
            continue
        end
        
        PHI_0(k,i) = U_inf*x;
        PSI_0(k,i) = U_inf*y;
        
        PHI_source(k,i) = q_source/(2*pi)*log(r1);
        PSI_source(k,i) = q_source/(2*pi)*th1;
        
        u_source(k,i) = q_source/(2*pi*r1)*cos(th1);
        v_source(k,i) = q_source/(2*pi*r1)*sin(th1);
        
        if  x > 0
            z1 = (D_net/2 + y)/(2*sqrt(epsilon_tau*x/U_inf)); 
            z2 = (D_net/2 - y)/(2*sqrt(epsilon_tau*x/U_inf));
            u1_pWake(k,i) = C_D_tot*U_inf/4*(erf(z1)+erf(z2));
        end
    end
    PHI = PHI_0+PHI_source;
    PSI = PSI_0+PSI_source;
    u = u_source+U_inf-u1_pWake;
    v = v_source;
end