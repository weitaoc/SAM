function [diff_u] = fun_diffusion_wus(cell_xyz,cell_R,Du,u,clv)

bio_parameters;

dist_cells = pdist2(cell_xyz,cell_xyz);
id_contact = (dist_cells<2*cell_R); % identify contacting cells according to the distance
id_contact = id_contact.*~eye(size(id_contact)); % remove the cell itself from contacting cells

%identify cells in layer 1
id_L1 = sqrt( cell_xyz(:,1).^2+cell_xyz(:,2).^2+cell_xyz(:,3).^2 )>=8.5;
id_L1c = repmat(id_L1,1,length(id_L1));
id_L1r = repmat(id_L1',length(id_L1),1);
id_L1mat = id_L1c.*id_L1r; % 1 for cell-cell pairs in L1

%identify cells in layer 2
id_L2 = ( sqrt( cell_xyz(:,1).^2+cell_xyz(:,2).^2+cell_xyz(:,3).^2 )>=7.5 ).*...
    ( sqrt( cell_xyz(:,1).^2+cell_xyz(:,2).^2+cell_xyz(:,3).^2 )<8.5 );
id_L2c = repmat(id_L2,1,length(id_L2));
id_L2r = repmat(id_L2',length(id_L2),1);
id_L2mat = id_L2c.*id_L2r;

%identify cells in deep layers
id_L3 = ones(size(cell_xyz,1),1);
id_L3 = id_L3-id_L1-id_L2;
id_L3c = repmat(id_L3,1,length(id_L3));
id_L3r = repmat(id_L3',length(id_L3),1);
id_L3mat = id_L3c.*id_L3r;

%construct clv over neighboring cells
clv_cmat = repmat(clv,1,length(clv));
clv_rmat = repmat(clv',length(clv),1);
clv_mat = 0.5 * ( clv_cmat + clv_rmat );
clv_mat = clv_mat .* id_contact .* (id_L1mat+id_L2mat+id_L3mat);

Lij_mat = dist_cells; % distance between contacting cells
Lij_mat = Lij_mat + diag(100*ones(size(u))); % avoid 0 in the distance
Aij_mat = pi * ( cell_R^2 - (Lij_mat/2).^2 ).*id_contact;
Uj_mat  = repmat(u',size(cell_xyz,1),1);
Ui_mat  = repmat(u,1,size(cell_xyz,1));
diff_u  = sum(Du./(1+(clv_mat/kcw3).^(4*n)).*Aij_mat.*(Uj_mat-Ui_mat)./Lij_mat,2);
