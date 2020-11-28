function b = wakeWidth(b0, C_Di, d_tw, x, xi, yi, Ni, N_elements, dy,UseFhSimCD)

b_incr = dy/100;
iter_max = 10;
b = b0;
N_cyl = Ni*N_elements;
I_wakePointCyl = find(xi(1:N_cyl/2)<x);
xi_m = xi(I_wakePointCyl);
yi_m = yi(I_wakePointCyl);
C_Di_m = C_Di(I_wakePointCyl);

for iter = 1:iter_max

    ui = sqrt(C_Di_m*d_tw./(x-xi_m)).*exp(-(b-yi_m).^2./(0.0888*C_Di_m*d_tw.*(x-xi_m)));
    u1_sCyl = 0.95*sum(ui);
    
    if u1_sCyl < 0.01
        return
    end
    b = b + b_incr;
end

end