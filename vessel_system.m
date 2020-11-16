function [eta_dot, ny_dot] = vessel_system(eta, ny, ny_r, Vc)
   % eta = [x, y, psi] ny = [u, v, r] Vc = [Vx, Vy, 0] 
   % coriolis_matrix = [C_RB(ny), C_A(ny_r)]
   % rotation_matrix = Rz(psi)
   
   %% Check if the constant matrices exists in Workspace
   ise = evalin( 'base', 'exist(''MRB'',''MA'', ''D'') == 1' );
   if ise
       disp('cool');
   end
       
   
end