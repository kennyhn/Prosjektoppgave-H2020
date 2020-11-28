function [u_sCyl, inWake, outWake] = CylModNetWake(x, y, xi, yi, C_Di, U, d_tw, lj, Ni, N_elements, cd, U2, D_net, b)

AsymptLimit = 100;
N_cyl = Ni*N_elements;
I_wakePointCyl = find(xi(1:N_cyl/2)<x);
xi_m = xi(I_wakePointCyl);
yi_m = yi(I_wakePointCyl);
C_Di_m = C_Di(I_wakePointCyl);
I_wakePointCyl = find((y-yi_m).^2 < AsymptLimit * 0.0888*C_Di_m.*d_tw.*(x-xi_m));
xi_m = xi_m(I_wakePointCyl);
yi_m = yi_m(I_wakePointCyl);
C_Di_m = C_Di_m(I_wakePointCyl);

% Sum of cylinders model
ui = sqrt(C_Di_m*d_tw./(x-xi_m)).*exp(-(y-yi_m).^2./(0.0888*C_Di_m*d_tw.*(x-xi_m)));
u1_sCyl = 0.95*U*sum(ui);
% Connecting to the potential flow by use of plane mixing layer theory
% Görtler (1942)

U1 = U-u1_sCyl;
%sigma = 13;%;*0.5/(U1/U2);
p1 = -1.473;
p2 = 2.6027;
% sigma = D_net.^(p1).*exp(p2).*D_net;
sigma = 13.5;

z1 = erf(sigma*((y-b)./D_net)/sqrt(x./D_net));
z2 = erf(sigma*((y+b)./D_net)/sqrt(x./D_net));

if ( z1 < 0.98 && z2 > -0.98) 
    inWake = 1;
    outWake = 0;
else
    inWake = 0;
    outWake = 1;
end
u_sCyl = U1 + (U2-U1)/2*(2+z1-z2);

end