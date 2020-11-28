function coefficientPlot(N_elements,C_N,C_D,C_T,C_L,phi,cl,cd)
    figure()
    plot(phi(1:N_elements/2)*180/pi, C_N(1:N_elements/2), phi(1:N_elements/2)*180/pi, C_T(1:N_elements/2))
    xlabel('$\phi$ [deg]','Interpreter', 'latex')
    legend('C_N', 'C_T')
    xlim([0 90])
    figure()
    plot(phi(1:N_elements/2)*180/pi, C_D(1:N_elements/2), phi(1:N_elements/2)*180/pi, C_L(1:N_elements/2))
    xlabel('$\phi$ [deg]','Interpreter', 'latex')
    legend('C_D', 'C_L')
    xlim([0 90])
    
    figure()
    hold on
    subplot(1,2,1)
    plot(phi(1:N_elements/2)*180/pi, C_D(1:N_elements/2)/cd)
    xlim([0 90])
    title('Normalized drag coefficient')
    xlabel('$\phi$ [deg]','Interpreter', 'latex')
    subplot(1,2,2)
    plot(phi(1:N_elements/2)*180/pi, C_L(1:N_elements/2)/cl)
    xlim([0 90])
    title('Normalized lift coefficient')
    xlabel('$\phi$ [deg]','Interpreter', 'latex')
    hold off

end