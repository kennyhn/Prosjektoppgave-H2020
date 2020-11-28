%**********************************************************************
%    Subroutine elementDiscretization
%    
%    Purpose: Discretizes a circular cylinder into straight 2D
%         line elements.
%    Method:  Uniform radial spacing 
%
%    Parameters
%    Input:
%    D_net      - Net diameter
%    N_elements - Number of elements
%    
%    Output:
%    x_nodes, y_nodes   - Vector of element coordinates
%    lj                 - Element length (equal for all)
%    phi                - Angle relative to element normal vector and x-axis
%
%    Programmed by: Lars Haug
%    Date: September 2020
%*********************************************************************

function [x_nodes, y_nodes, lj, theta, phi] = elementDiscretization(N_elements, D_net)
N_nodes = N_elements+1;
dtheta = 2*pi/N_elements;
theta = [0:dtheta:N_elements*dtheta]-N_elements/2*dtheta/2;
[x_nodes, y_nodes] = pol2cart(theta, ones(1,N_nodes)*D_net/2);
% x_nodes = zeros(1,N_elements+1);
% y_nodes = -D_net/2:D_net/(N_elements):D_net/2;

lj = sqrt((x_nodes(end)-x_nodes(end-1))^2 + (y_nodes(end)-y_nodes(end-1))^2);

phi = zeros(1, N_elements);

for e = 1:N_elements
    elem_centre_x = (x_nodes(e+1)+x_nodes(e))/2;
    elem_centre_y = (y_nodes(e+1)+y_nodes(e))/2;
    
    phi(e) = cart2pol(elem_centre_x,elem_centre_y);
    
end
end