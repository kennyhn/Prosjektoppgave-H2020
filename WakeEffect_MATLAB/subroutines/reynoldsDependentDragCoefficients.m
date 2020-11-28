function [C_N, C_T, C_D, C_L, cd, cl] = reynoldsDependentDragCoefficients(Re, Sn, phi) 
    
    phi = phi - pi*(phi<pi & phi>pi/2) + pi*(phi>-pi & phi<-pi/2);
    
    if Re < 10^(3/2)
        C_D_cirCyl = 1 + 10/(Re^(2/3));
    else
        C_D_cirCyl = -78.46675 + 254.73873*log10(Re) - 327.8864*(log10(Re))^2 + 223.64557*(log10(Re))^3 - 87.92234*(log10(Re))^4 ...
            + 20.00769*(log10(Re))^5 - 2.44894*(log10(Re))^6 + 0.12497*(log10(Re))^7;
    end
    
    C_D_Loland = 0.04 + (-0.04 + 0.33*Sn + 6.54*Sn^2 - 4.88*Sn^3)*cos(phi);
    
    % Fourier Components  
    C_N_45 = C_D_cirCyl*Sn*(2-Sn)/(2*(1-Sn)^2)*(cos(pi/4))^2;
    C_T_45 = pi/4*4*C_N_45/(8+C_N_45);  
    cd = C_D_cirCyl*Sn*(2-Sn)/(2*(1-Sn)^2);
    cl = (0.5*cd-C_T_45)/sqrt(2);
    a3 = 0.1;
    a1 = 1-a3;
    b4 = 0.1;
    b2 = 1;

    C_N = C_D_cirCyl*Sn*(2-Sn)/(2*(1-Sn)^2)*(cos(phi))^2;
    C_T = phi*4*C_N/(8+C_N);
    C_D = cd*(a1*cos(phi) + a3*cos(3*phi));
    C_L = cl*(b2*sin(2*phi) + b4*sin(4*phi));
    %C_D = C_D_Loland;
end