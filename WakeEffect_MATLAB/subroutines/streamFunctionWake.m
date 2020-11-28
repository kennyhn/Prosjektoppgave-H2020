function streamFunctionWake(x_nodes, y_nodes,  theta_nodes, N_elements, x_vec, y_vec)
%{
    *************************************************************
    ARGUMENTS
    x_nodes:        horizontal coordinate of nodes
    y_nodes:        vertical coordinate of nodes
    theta_nodes:    angle to node from origin (0-2pi)
    ---------
    VARIABLES
    g              horizontal coordinate of element j midpoint
    gdot           slope of element j as a function of eta
    eta_j          vertical coordinate of element j midpoint
    deta_j         vertical length component of element j as a function of eta
    nStart         node start index of curve section
    nEnd           node end index of curve section
    kappa          pressure coefficient for flow through screen
    alfa_j         attack angle of element midpoint
    ---------
    *************************************************************
%}

N_nodes = N_elements+1;
nStart = (N_nodes+1)/2;
nEnd = N_nodes;%-(nStart-1)/2;
%nStart+1;

% Upstream screen
alfa_j = pi-(theta_nodes(nStart:nEnd-1)+theta_nodes(nStart+1:nEnd))/2;
eta_j = (y_nodes(nStart:nEnd-1)+y_nodes(nStart+1:nEnd))/2;
deta_j = y_nodes(nStart:nEnd-1)-y_nodes(nStart+1:nEnd);
g = (x_nodes(nStart:nEnd-1)+x_nodes(nStart+1:nEnd))/2;
gdot = (x_nodes(nStart:nEnd-1)-x_nodes(nStart+1:nEnd))./deta_j;

qmax = 0.3/(N_elements/2);
qmin = 0.0001;
q = [qmin:(qmax-qmin)/(N_elements/4-1):qmax qmax:-(qmax-qmin)/(N_elements/4-1):qmin];
%q = ones(1,N_elements/2)*qmax;
iterMax = 10;
tolerance = 1e-2;
nel = N_elements/2;

%for iter = 1:iterMax
for j = 1:nel
    z = (g(j)+eps)+eta_j(j)*1i;
    ln1 = -1./(gdot+1i).*log(sin(pi*(z-(g+1i*eta_j(j)))/(2i)));
    ln2 = -1./(gdot-1i).*log(sin(pi*(z-(g+1i*eta_j(j)))/(2i)));
    uv1(j) = 1 + sum( q/(2*pi).*sqrt(gdot.^2+1).*(ln1 + ln2) );
end
u1 = real(uv1);
v1 = imag(uv1);
u2 = 1 + 0.5*sum(q.*sqrt(gdot.^2+1).*deta_j);

UNM = u1.*cos(alfa_j+v1.*sin(alfa_j));
%UNP = E.*UNM
%end

for n=1:length(x_vec)
    zeta = x_vec(n);
    for m=1:length(y_vec)
        eta = y_vec(m);
        z = zeta + eta*1i;
        
        ln1 = -1./(gdot+1i).*log(sin(pi*(z-(g+1i*eta_j))/(2i)));
        ln2 = -1./(gdot-1i).*log(sin(pi*(z-(g-1i*eta_j))/(2i)));
        uv1(m,n) = 0 + sum( q/(2*pi).*sqrt(gdot.^2+1).*(ln1 + ln2) );
        
        ln1 = log(sin(pi*(z-(g+1i*eta_j))/(2i)));
        ln2 = log(sin(pi*(z-(g-1i*eta_j))/(2i)));
        W(m,n) = 0 + sum( q/(2*pi).*sqrt(gdot.^2+1).*(ln1 + ln2) ); 
    end
end
u = real(uv1);
v = imag(uv1);
PHI = real(W);
PSI = imag(W);
[X,Y] = meshgrid(x_vec, y_vec);

figure()
hold on
plot(g,eta_j,'*')
plot(x_nodes(nStart:nEnd), y_nodes(nStart:nEnd),'-')
scale = 0.5;
quiver(X,Y,u,v,scale)
contour(X,Y,PSI);
hold off
end