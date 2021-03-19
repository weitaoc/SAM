function [diff_u] = fun_diffusion(cell_xyz,cell_R,Du,u)

dist_cells = pdist2(cell_xyz,cell_xyz);
id_contact = (dist_cells<2*cell_R); % identify contacting cells according to the distance
id_contact = id_contact.*~eye(size(id_contact)); % remove the cell itself from contacting cells

Lij_mat = dist_cells; % distance between contacting cells
Lij_mat = Lij_mat + diag(100*ones(size(u))); % avoid 0 in the distance
Aij_mat = pi * ( cell_R^2 - (Lij_mat/2).^2 ).*id_contact;
Uj_mat  = repmat(u',size(cell_xyz,1),1);
Ui_mat  = repmat(u,1,size(cell_xyz,1));
diff_u  = sum(Du*Aij_mat.*(Uj_mat-Ui_mat)./Lij_mat,2);
