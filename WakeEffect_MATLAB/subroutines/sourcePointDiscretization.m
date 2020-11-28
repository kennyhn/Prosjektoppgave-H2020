%**********************************************************************
%    Subroutine sourcePointDiscretization
%    
%    Purpose: Distributes source points along an element to represent
%             a twine.
%    Method:  Uniform radial spacing 
%
%    Parameters
%    Input:
%    x_nodes, y_nodes   - Vector of element coordinates
%    N_elements         - Number of elements
%    Ni                 - Number of source points in single element 
%    
%    Output:
%    xi, yi             - Vector of source point coordinates
%
%    Programmed by: Lars Haug
%    Date: September 2020
%*********************************************************************

function [xi, yi, exi, eyi, phi] = sourcePointDiscretization(x_nodes, y_nodes, N_elements, Ni)
% Memory Allocation 
phi = zeros(1, N_elements);
yi = zeros(1, N_elements*Ni);
xi = zeros(1, N_elements*Ni);

for e = 1:N_elements
    lambda_x = (x_nodes(e+1)-x_nodes(e))/Ni;
    lambda_y = (y_nodes(e+1)-y_nodes(e))/Ni;

    exi = (x_nodes(e+1)+x_nodes(e))/2;
    eyi = (y_nodes(e+1)+y_nodes(e))/2;

    phi(e) = cart2pol(exi,eyi);


    if lambda_y == 0
        yi(1,(e-1)*Ni+1:e*Ni) = y_nodes(e).*ones(1,Ni);
    else
        yi(1,(e-1)*Ni+1:e*Ni) = y_nodes(e)+lambda_y/2:lambda_y:(y_nodes(e+1)-lambda_y/2);
    end

    if lambda_x == 0
        xi(1,(e-1)*Ni+1:e*Ni) = x_nodes(e).*ones(1,Ni);
    else
        xi(1,(e-1)*Ni+1:e*Ni) = x_nodes(e)+lambda_x/2:lambda_x:(x_nodes(e+1)-lambda_x/2);
    end
end

end