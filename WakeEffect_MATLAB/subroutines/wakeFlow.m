%**********************************************************************
%    Subroutine wakeFlow
%    
%    Purpose: Calculation of velocities in flow field
%      
%    Method:  Mean turbulent wake velocites are calculated using 
%         linearized turbulent wake equations at individual cylinders 
%         that resemples the twines. Velocities in front of cage and  
%         outside of cage are calculated using lagally's theorem.  
%    Parameters
%    Input:
%    rho_s          - Water density
%    U_dnet         - Inflow velocity of downstream net sections
%    U_inf          - Uniform ambient velocity
%    Sn             - Solidity ratio
%    D_net          - Net diameter
%    d_tw           - Twine diameter
%    N_elements     - Number of screen elements  
%    C_D            - Vector of drag coefficients on screen elements
%    x_vec, y_vec   - Elements coordinate 
%    xi, yi         - Source points coordinate
%
%    Output:
%    u, v           - Velocitiessss
%
%    Programmed by: Lars Haug
%    Date: September 2020
%*********************************************************************

function [u, v, PHI, PSI, u1_pWake] = wakeFlow(U_dnet, D_net, N_elements, U_inf, rho_s, C_D, lj, x_vec, y_vec, xi, yi, C_Di, d_tw, Ni, cd, dy, ymax, method)

error_param = 1;

% Memory allocation
PHI_0 = zeros(length(x_vec),length(y_vec));
PSI_0 = zeros(length(x_vec),length(y_vec));

PHI_source = zeros(length(x_vec),length(y_vec));
PSI_source = zeros(length(x_vec),length(y_vec));

u_source = zeros(length(x_vec),length(y_vec));
v_source = zeros(length(x_vec),length(y_vec));

u1_pWake = zeros(length(x_vec),length(y_vec));
PSI_wake = zeros(length(x_vec),length(y_vec));
PHI_wake = zeros(length(x_vec),length(y_vec));

u1_sCyl = zeros(length(x_vec),length(y_vec));
u_sCyl = zeros(length(x_vec),length(y_vec));
inWake = zeros(length(x_vec),length(y_vec));
outWake = ones(length(x_vec),length(y_vec));
b = D_net/2;

U_elem = [ones(1,N_elements/2)*U_dnet ones(1,N_elements/2)*U_inf];
F_D = 0.5*rho_s.*U_elem.^2.*C_D*lj;

C_D_tot = 2*sum(F_D)/(rho_s*U_inf^2*D_net);
epsilon_tau = error_param*0.0222*C_D_tot*D_net*U_inf;

q_source = sum(F_D)/(rho_s*U_inf);

% Loading variables
N_operations = length(x_vec);
count = -1;
indx = 1;
LoadCount = N_operations*[0.1:0.05:1];
for k = 1:length(x_vec)
    x = x_vec(k);
    if x > 0 && method == 1
        b = wakeWidth(b, C_Di, d_tw, x, xi, yi, Ni, N_elements, dy);
        [th1, r1] = cart2pol(x,b);
        u_source_boundary = U_inf + q_source/(2*pi*r1)*cos(th1);
    end
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
        
        u_source(k,i) = U_inf + q_source/(2*pi*r1)*cos(th1);
        v_source(k,i) = q_source/(2*pi*r1)*sin(th1);
        if x > 0
            if method == 1
                [u_sCyl(k,i), inWake(k,i), outWake(k,i)] = CylModNetWake(x, y, xi, yi, C_Di, U_dnet, d_tw, lj, Ni, N_elements, cd, u_source_boundary, D_net, b);
                
                PSI_wake(k,i) =  u_sCyl(k,i)*y;
                PHI_wake(k,i) =  u_sCyl(k,i)*x;
                u1_sCyl(k,i) = U_inf - u_sCyl(k,i);
                
            elseif method == 2
                z1 = (D_net/2 + y)/(2*sqrt(epsilon_tau*x/U_inf));
                z2 = (D_net/2 - y)/(2*sqrt(epsilon_tau*x/U_inf));
                if ( erf(z1) + erf(z2) > 0.02 ) 
                    inWake(k,i) = 1;
                    outWake(k,i) = 0;
                else
                    inWake(k,i) = 0;
                    outWake(k,i) = 1;
                end
                u1_pWake(k,i) = C_D_tot*U_inf/4*(erf(z1)+erf(z2));
                PSI_wake(k,i) = -u1_pWake(k,i)*y;
                PHI_wake(k,i) = -u1_pWake(k,i)*x;
            else
                disp("Error: Choose a method")
            end
        end
        
    end
    % Updating Loading Screen
    count = count + 1;
    if count == int64(LoadCount(indx))
        disp(['Progress: ' num2str(5*indx) ' %'])
        indx = indx+1;
    end
end


if method == 1
    PHI = PHI_0+PHI_source;
    PSI = (PSI_0+PSI_source).*outWake+PSI_wake.*inWake;
    u = u_source.*outWake + u_sCyl.*inWake;
    v = v_source.*outWake;
    
elseif method == 2
    PHI = PHI_0+PHI_source;
    PSI = PSI_0+PSI_source+PSI_wake;
    u = u_source.*outWake+(U_inf-u1_pWake).*inWake;
    v = v_source.*outWake;
end

disp(['Progress: 100 %'])
end