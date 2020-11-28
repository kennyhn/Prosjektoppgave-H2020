function flowPlot(xmin, xmax, ymax, D_net, dx, dy, x_vec, y_vec, x_nodes, y_nodes, U_dnet, U_inf, u, v, PSI, cd, cd_loland, PaperComparisson, boundary, savefigure)

% delta_1 = zeros(1,length(x_vec));
% 
% for k = 201:length(x_vec)
%     delta_1(k) = (trapz(u1_pWake(k,:)./U_inf))/D_net;
% end
u1_loland = (1-0.46*cd_loland)^2*ones(1,length(y_vec));

obstacle = zeros(length(x_vec),length(y_vec));
[X,Y] = meshgrid(x_vec,y_vec);

kmin = ceil(abs(xmin+D_net/2)/dx)+1;
kmax = floor(abs(xmin-D_net/2)/dx)+1;
imin = ceil(abs(ymax-D_net/2)/dy)+1;
imax = floor(abs(ymax+D_net/2)/dy)+1;
for k=kmin:kmax
    x = x_vec(k);
    for i=imin:imax
        y = y_vec(i);
        if x^2+y^2 < (D_net/2)^2
            obstacle(k,i) = 1;
            continue
        end
        
    end
end

u_R = (u+obstacle.*(U_dnet))/U_inf;
v_R = v/U_inf;

% figure()
% plot(x_vec./D_net, delta_1)
% xlabel('x-coordinate')
% ylabel('$\delta_1$','Interpreter','latex')
% xlim([1 5])
% title('Discplacement Thickness')


figure()
hold on
dx1 = 0.005; dx2 = 0.05; dx3 = 0.01;
[~,c] = contourf(X./D_net,Y./D_net,u_R',[0.5:dx1:U_dnet/U_inf*0.9 U_dnet/U_inf*0.9+dx2:dx2:U_dnet/U_inf*1.25 U_dnet/U_inf*1.25+dx3:dx3:1.2]);
pbaspect([xmax-xmin 2*ymax 1])
colormap(jet);
c.LineWidth = 0.0001;
c.LineColor = [0.35 0.35 0.35];
plot(x_nodes./D_net,y_nodes./D_net,'r','LineWidth',3)
%plot(boundary(:,1)./D_net,boundary(:,2)./D_net,'--w','LineWidth',3)
xlabel('x-axis, [$\frac{x}{D_{net}}$]','Interpreter','latex')
ylabel('y-axis, [$\frac{y}{D_{net}}$]','Interpreter','latex')
colorbar
[~,~,u_axis] = find(u_R);
caxis([min(u1_loland) max(u_axis)])
clear u_axis
legend('Reduced Velocity, $\frac{u}{U_{\infty}}$','Circular Cage','Interpreter','latex')
title({'Horizontal Velocity Distribution for the Sum of Cylinder Model',['D_{net} = ' num2str(D_net) ', U_{\infty} = ' num2str(U_inf)]})
if savefigure == true
    savefig('fig\Contour_u.fig')
    saveas(gcf,'fig\Contour_u.png')
end
hold off

figure()
hold on
dx1 = 0.005; dx2 = 0.001; dx3 = 0.005; wakelim = 0.001;
[~,c] = contourf(X./D_net,Y./D_net,v_R',[-0.15:dx1:-wakelim -wakelim+dx2:dx2:wakelim wakelim+dx3:dx3:0.15]);
pbaspect([xmax-xmin 2*ymax 1])
colormap(jet);
c.LineWidth = 0.01;
c.LineColor = [0.35 0.35 0.35];
plot(x_nodes./D_net,y_nodes./D_net,'r','LineWidth',3)
xlabel('x-axis, [$\frac{x}{D_{net}}$]','Interpreter','latex')
ylabel('y-axis, [$\frac{y}{D_{net}}$]','Interpreter','latex')
colorbar
caxis([min(min(v_R)) max(max(v_R))])
legend('Reduced vertical Velocity, $\frac{v}{U_{\infty}}$','Circular Cage','Interpreter','latex')
title({'Vertical Velocity Distribution for the Sum of Cylinder Model',['D_{net} = ' num2str(D_net) ', U_{\infty} = ' num2str(U_inf)]})
if savefigure == true    
    savefig('fig\Contour_v.fig')
    saveas(gcf,'fig\Contour_v.png')
end
hold off

U_flux = sqrt((v_R*U_inf).^2+(u_R*U_inf).^2)';
Flux1 = 1025*trapz(Y(:,1),U_flux(:,1));
Flux2 = 1025*trapz(Y(:,end),U_flux(:,end));
% Flux1-Flux2
% (Flux1-Flux2)/Flux1
figure()
hold on
scale = 0.5;
XSVc = int64(1:20:(xmax-xmin)/dx);
YSVc = int64(1:20:(2*ymax)/dx);
[~, c] = contour(X./D_net,Y./D_net,PSI',32);
c.LineWidth = 1.7;
colormap(gray)
quiver(X(YSVc,XSVc)./D_net,Y(YSVc,XSVc)./D_net,u_R(XSVc,YSVc)',v_R(XSVc,YSVc)',scale)
pbaspect([xmax-xmin 2*ymax 1])
plot(x_nodes./D_net,y_nodes./D_net,'r','LineWidth',3)
xlabel('x-axis, [$\frac{x}{D_{net}}$]','Interpreter','latex')
ylabel('y-axis, [$\frac{y}{D_{net}}$]','Interpreter','latex')

legend('Streamlines','Velocity Vectors')
hold off


xMeasure = [-2 -1.5 -1.0 1 1.5 3]*D_net;
IMeasure = round(-xmin/dx+xMeasure/dx+1,0);

figure()
hold on
lgd = {};
color = [255,128,0; 102,0,204; 102,255,255; 204,0,0; 0,153,0; 0,0,204]./255;
for i = 1:length(IMeasure)
    if xMeasure(i)<0
        style = '--';
    else
        style = '-';
    end
    plot(y_vec./D_net,u_R(IMeasure(i),:),style, 'Color' ,color(i,:))
    lgd{1,i} = ['x/D_{net}=' num2str(xMeasure(i)/D_net)];
end


if PaperComparisson == 1
    u_paper1 = readContourPlot("data\turner_etal_velocontour_exit.png",[-0.91 0.91],[1.39 0.16],[0.1864 1.2725],false);
    [u_paper2, pixels_exit] = readContourPlot("data\turner_etal_velocontour_1d.png",[-0.91 0.91],[1.39 0.16],[0.1864 1.2725],false);
    u_paper3 = readContourPlot("data\turner_etal_velocontour_2d.png",[-0.91 0.91],[1.39 0.16],[0.1864 1.2725],false); 
    incr = (1.39-0.16)/pixels_exit(1);
    row = round((0.5-0.16)/incr);
    xvals = -0.91:(2*0.91)/(pixels_exit(2)-1):0.91;

    plot(xvals,u_paper1(row,1:pixels_exit(2)),'.')
    plot(xvals,u_paper2(row,1:pixels_exit(2)),'.')
    plot(xvals,u_paper3(row,1:pixels_exit(2)),'.')
    lgd{1,length(IMeasure)+1} = 'Turner et. al, x/D_{net}= 0.75';
    lgd{1,length(IMeasure)+2} = 'Turner et. al, x/D_{net}= 1.5';
    lgd{1,length(IMeasure)+3} = 'Turner et. al, x/D_{net}= 2.5';
    i = i + 3;
end

if PaperComparisson == 2
    T = readmatrix('data\gansel_etal_u_0_2_Sn_25_yD_1');
    u_paper1 = (0.2-T(:,2))./0.2; xvals1 = T(:,1);
    T = readmatrix('data\gansel_etal_u_0_2_Sn_25_yD_1_5');
    u_paper2 = (0.2-T(:,2))./0.2; xvals2 = T(:,1);
    T = readmatrix('data\gansel_etal_u_0_2_Sn_25_yD_3');
    u_paper3 = (0.2-T(:,2))./0.2; xvals3 = T(:,1);

    plot(xvals1,u_paper1,'*', 'Color', color(4,:))
    plot(xvals2,u_paper2,'d', 'Color', color(5,:))
    plot(xvals3,u_paper3,'o', 'Color', color(6,:))
    lgd{1,length(IMeasure)+1} = 'Gansel et. al, x/D_{net}= 1';
    lgd{1,length(IMeasure)+2} = 'Gansel et. al, x/D_{net}= 1.5';
    lgd{1,length(IMeasure)+3} = 'Gansel et. al, x/D_{net}= 3';
    i = i + 3;
end
    plot(y_vec./D_net,u1_loland,'--k')
    lgd{1,i+1} = 'Lølands formula: (1-0.46cd)^2';

    legend(lgd,'Location','southeast')
    
    xlabel('y-axis, [$\frac{y}{D_{net}}$]','Interpreter','latex')
    ylabel('[$\frac{u}{U_{\infty}}$]','Interpreter','latex')
    xlim([0 ymax/D_net])
    title({'Transverse Velocity Profile at specified distances',['D_{net} = ' num2str(D_net) ', U_{\infty} = ' num2str(U_inf)]})
    if savefigure == true
        savefig('fig\Vel_profile_u.fig')
        saveas(gcf,'fig\Vel_profile_u.png')
    end
    hold off

end

