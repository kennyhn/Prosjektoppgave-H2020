function [u1_mCyl, u1_sCyl] = CylModNetWake(x, y, xi, yi, C_Di, U, d_tw, L, Ni, N_elements, cd)



N_cyl = Ni*N_elements;
I_wakePointCyl = find(xi(1:N_cyl/2)<x);
xi_m = xi(I_wakePointCyl);
yi_m = yi(I_wakePointCyl);
C_Di_m = C_Di(I_wakePointCyl);

% Sum of cylinders model
ui = sqrt(C_Di_m*d_tw./(x-xi_m)).*exp(-(y-yi_m).^2./(0.0888*C_Di_m*d_tw.*(x-xi_m)));
u1_sCyl = 0.95*U*sum(ui)*(0.46*cd);

% Modified sum of cylinders model
if isempty(xi_m)
    u1_mCyl = 0;
else
    C_D_thr = 1.2;
    epsilon_zero = 0.0222.*C_D_thr.*d_tw.*U;
    
    a = epsilon_zero.*(x-xi_m).*(1+(Ni-1).*tanh((x-xi_m)./C_Di_m*d_tw));
    A = sqrt(pi*U./a);
    B = exp(-(y-yi_m).^2.*U./(4*a));
    ui = C_Di_m.*L.*U.*A.*B;
    u1_mCyl = 1/(4*pi).*U*sum(ui);
end

end