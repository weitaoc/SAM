% this code compute WUS-CLV3-Cytokinin system in 3D: the 3D tissue template
% is borrowed from 2013 MSB paper with 1366 spherical cells at radius 1.
% Overlapping cells are considered to contact with each other and diffusing
% molecules will follow passive transport between contacting cells
clear
close all
global EL WB
EL = 1;
WB = 1;

bio_parameters;

load WCK.init
cell_R = 1;
cell_xyz = WCK(:,1:3);
Num_cell = size(cell_xyz,1);


% index contacting cells
dist_cells = pdist2(cell_xyz,cell_xyz);
id_contact = (dist_cells<2*cell_R); % identify contacting cells according to the distance
id_contact = id_contact.*~eye(size(id_contact)); % remove the cell itself from contacting cells
figure(1); hold on; % show overlapping cells for example
temp = 1000;
scatter3(cell_xyz(temp,1),cell_xyz(temp,2),cell_xyz(temp,3),1000,'o','filled')
scatter3(cell_xyz(id_contact(temp,:)>0,1),cell_xyz(id_contact(temp,:)>0,2),cell_xyz(id_contact(temp,:)>0,3),1000,'ro')
hold off; 

% index cells at bottom boundary: choose a circular plane below the tissue,
% randomly select large number of points within the plane, bottom of the 
% tissue is defined to be the cells which give the shortest distance from 
% a point on the plane to the tissue
zmin = min(cell_xyz(:,3));
r_bottom = 10;
Nrand = 10000;
x_rand_bottom = r_bottom*(2*rand(Nrand,1)-1);
y_rand_bottom = r_bottom*(2*rand(Nrand,1)-1);
ind = find(x_rand_bottom.^2+y_rand_bottom.^2<r_bottom^2);
x_rand_bottom = x_rand_bottom(ind);
y_rand_bottom = y_rand_bottom(ind);
z_rand_bottom = (zmin-2)*ones(size(x_rand_bottom));
% scatter(x_rand_bottom,y_rand_bottom,'o');
dist_plane_cells = pdist2([x_rand_bottom, y_rand_bottom, z_rand_bottom],cell_xyz);
[val,ind] = min(dist_plane_cells,[],2);
id_bottom = unique(ind);
figure(2); hold on; %identify cells at bottom
id_above = 1:Num_cell;
id_above(id_bottom) = [];
scatter3(cell_xyz(id_above,1),cell_xyz(id_above,2),cell_xyz(id_above,3),500,'o','filled')
scatter3(cell_xyz(id_bottom,1),cell_xyz(id_bottom,2),cell_xyz(id_bottom,3),500,'ro')
hold off;

% define variables
wusR = zeros(Num_cell,1);
wn   = wusR;
wc   = wusR;
clv  = wusR;
ck   = wusR;
ckR  = wusR;
CK   = wusR;

% define source function
% wusR_source = @(x,y,z) A1*exp(-(x-xc1).^2/sigmax1.^2-(y-yc1).^2/sigmay1.^2-(z-zc1).^2/sigmaz1.^2).*(x.^2+y.^2+z.^2<=L_ckR^2);
wusR_source = @(x,y,z) A1*(z<=L_wus).*(x.^2+y.^2<=L_clv^2);
% clv_source  = @(x,y,z) A2*exp(-(x-xc2).^2/sigmax2.^2-(y-yc2).^2/sigmay2.^2-(z-zc2).^2/sigmaz2.^2).*(x.^2+y.^2+z.^2>L_ckR^2);
clv_source  = @(x,y,z) A2.*(x.^2+y.^2<=L_clv^2);
ck_source   = @(x,y,z) A3*exp(-(x-xc3).^2/sigmax3.^2-(y-yc3).^2/sigmay3.^2-(z-zc3).^2/sigmaz3.^2).*(x.^2+y.^2+z.^2<=L_ckR^2);
ckR_source  = @(x,y,z) A4*(x.^2+y.^2+z.^2<=L_ckR^2);
figure(3); % show expression domain
% cross section identification
id_cs = find(abs(cell_xyz(:,1))<1);
subplot(4,2,1)
scatter3(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3),50,wusR_source(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3)),'filled');
title('wus mRNA');
subplot(4,2,2)
scatter3(cell_xyz(id_cs,1),cell_xyz(id_cs,2),cell_xyz(id_cs,3),50,wusR_source(cell_xyz(id_cs,1),cell_xyz(id_cs,2),cell_xyz(id_cs,3)),'filled');
view([90 0]); colorbar
subplot(4,2,3)
scatter3(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3),50,clv_source(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3)),'filled');
title('CLV3');
subplot(4,2,4)
scatter3(cell_xyz(id_cs,1),cell_xyz(id_cs,2),cell_xyz(id_cs,3),50,clv_source(cell_xyz(id_cs,1),cell_xyz(id_cs,2),cell_xyz(id_cs,3)),'filled');
view([90 0]);colorbar
subplot(4,2,5)
scatter3(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3),50,ck_source(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3)),'filled');
title('Cytokinin');
subplot(4,2,6)
scatter3(cell_xyz(id_cs,1),cell_xyz(id_cs,2),cell_xyz(id_cs,3),50,ck_source(cell_xyz(id_cs,1),cell_xyz(id_cs,2),cell_xyz(id_cs,3)),'filled');
view([90 0]);colorbar
subplot(4,2,7)
scatter3(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3),50,ckR_source(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3)),'filled');
title('Cytokinin receptor');
subplot(4,2,8)
scatter3(cell_xyz(id_cs,1),cell_xyz(id_cs,2),cell_xyz(id_cs,3),50,ckR_source(cell_xyz(id_cs,1),cell_xyz(id_cs,2),cell_xyz(id_cs,3)),'filled');
view([90 0]);colorbar


T = 300;
Tplot = 0.1;
Tsave = 30;
dt = 0.01;
Nt = floor(T/dt);
Np = floor(Tplot/dt);
Ns = floor(Tsave/dt);


for it = 1:Nt
    
    % update by RK2
    
    [diff_wc] = fun_diffusion_wus(cell_xyz,cell_R,Dw,wc,clv);
    [diff_clv] = fun_diffusion(cell_xyz,cell_R,Dc,clv);
    [diff_ck] = fun_diffusion(cell_xyz,cell_R,Dck,ck);
    [fwusR,fwn,fwc,fclv,fck,fckR,fCK] = fun_reaction(cell_xyz,wusR,wn,wc,clv,ck,ckR,CK,wusR_source,clv_source,ck_source,ckR_source);
    wusR_temp = wusR + dt * fwusR;
    wn_temp   = wn   + dt * fwn;
    wc_temp   = wc   + dt * ( diff_wc + fwc );
    clv_temp  = clv  + dt * ( diff_clv + fclv );
    ck_temp   = ck   + dt * ( diff_ck + fck );
    ckR_temp  = ckR  + dt * fckR;
    CK_temp   = CK   + dt * fCK;
    %impose leaky boundary condition at bottom
    wc_temp(id_bottom)  = wc_temp(id_bottom)  - falpha * wc_temp(id_bottom);
    clv_temp(id_bottom) = clv_temp(id_bottom) - falpha * clv_temp(id_bottom);
    ck_temp(id_bottom)  = ck_temp(id_bottom)  - falpha * ck_temp(id_bottom);
    
    
    [diff_wc] = fun_diffusion_wus(cell_xyz,cell_R,Dw,wc_temp,clv_temp);
    [diff_clv] = fun_diffusion(cell_xyz,cell_R,Dc,clv_temp);
    [diff_ck] = fun_diffusion(cell_xyz,cell_R,Dck,ck_temp);
    [fwusR,fwn,fwc,fclv,fck,fckR,fCK] = fun_reaction(cell_xyz,wusR_temp,wn_temp,wc_temp,clv_temp,ck_temp,ckR_temp,CK_temp,wusR_source,clv_source,ck_source,ckR_source);
    wusR = 0.5 * wusR + 0.5 * wusR_temp + 0.5 * dt * fwusR;
    wn   = 0.5 * wn   + 0.5 * wn_temp   + 0.5 * dt * fwn;
    wc   = 0.5 * wc   + 0.5 * wc_temp   + 0.5 * dt * ( diff_wc + fwc );
    clv  = 0.5 * clv  + 0.5 * clv_temp  + 0.5 * dt * ( diff_clv + fclv );
    ck   = 0.5 * ck   + 0.5 * ck_temp   + 0.5 * dt * ( diff_ck + fck );
    ckR  = 0.5 * ckR  + 0.5 * ckR_temp  + 0.5 * dt * fckR;
    CK   = 0.5 * CK   + 0.5 * CK_temp   + 0.5 * dt * fCK;
    %impose leaky boundary condition at bottom
    wc(id_bottom)  = wc(id_bottom)  - falpha * wc(id_bottom);
    clv(id_bottom) = clv(id_bottom) - falpha * clv(id_bottom);
    ck(id_bottom)  = ck(id_bottom)  - falpha * ck(id_bottom);
    
    % nonnegativity
    wusR = wusR.*(wusR>0);
    wn   = wn.*(wn>0);
    wc   = wc.*(wc>0);
    clv  = clv.*(clv>0);
    ck   = ck.*(ck>0);
    ckR  = ckR.*(ckR>0);
    CK   = CK.*(CK>0);
    
    if mod(it,Np) == 0
        figure(4);
        id_cs = find(abs(cell_xyz(:,1))<1);
        subplot(7,2,1)
        scatter3(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3),50,wusR,'filled');
        title('wus mRNA');
        subplot(7,2,2)
        scatter3(cell_xyz(id_cs,1),cell_xyz(id_cs,2),cell_xyz(id_cs,3),50,wusR(id_cs),'filled');
        view([90 0]); colorbar
        subplot(7,2,3)
        scatter3(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3),50,wn,'filled');
        title('wus in nuclei');
        subplot(7,2,4)
        scatter3(cell_xyz(id_cs,1),cell_xyz(id_cs,2),cell_xyz(id_cs,3),50,wn(id_cs),'filled');
        view([90 0]); colorbar
        subplot(7,2,5)
        scatter3(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3),50,wc,'filled');
        title('wus in cytoplasm');
        subplot(7,2,6)
        scatter3(cell_xyz(id_cs,1),cell_xyz(id_cs,2),cell_xyz(id_cs,3),50,wc(id_cs),'filled');
        view([90 0]); colorbar
        subplot(7,2,7)
        scatter3(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3),50,clv,'filled');
        title('CLV3');
        subplot(7,2,8)
        scatter3(cell_xyz(id_cs,1),cell_xyz(id_cs,2),cell_xyz(id_cs,3),50,clv(id_cs),'filled');
        view([90 0]);colorbar
        subplot(7,2,9)
        scatter3(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3),50,ck,'filled');
        title('Cytokinin');
        subplot(7,2,10)
        scatter3(cell_xyz(id_cs,1),cell_xyz(id_cs,2),cell_xyz(id_cs,3),50,ck(id_cs),'filled');
        view([90 0]);colorbar
        subplot(7,2,11)
        scatter3(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3),50,ckR,'filled');
        title('Cytokinin receptor');
        subplot(7,2,12)
        scatter3(cell_xyz(id_cs,1),cell_xyz(id_cs,2),cell_xyz(id_cs,3),50,ckR(id_cs),'filled');
        view([90 0]);colorbar
        subplot(7,2,13)
        scatter3(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3),50,CK,'filled');
        title('Cytokinin');
        subplot(7,2,14)
        scatter3(cell_xyz(id_cs,1),cell_xyz(id_cs,2),cell_xyz(id_cs,3),50,CK(id_cs),'filled');
        view([90 0]);colorbar
        drawnow
    end
    
    if mod(it,Ns) == 0
        save(['data_T' num2str(it*dt) '.mat']);
        fid = fopen(['wus_data3d_T' num2str(it*dt) '.vtk'],'w+');
        fprintf(fid,'%s\n','# vtk DataFile Version 3.0');
        fprintf(fid,'%s\n','Points representing individual cells');
        fprintf(fid,'%s\n\n','ASCII');
        fprintf(fid,'%s\n','DATASET UNSTRUCTURED_GRID');
        fprintf(fid,'%s ','POINTS');
        fprintf(fid,'%d ',Num_cell);
        fprintf(fid,'%s\n','float');
        for i = 1:Num_cell
            fprintf(fid,'%f %f %f\n',cell_xyz(i,1),cell_xyz(i,2),cell_xyz(i,3));
        end
        
        fprintf(fid,'\n');
        fprintf(fid,'%s ','CELLS');
        fprintf(fid,'%d ',Num_cell);
        fprintf(fid,'%d\n',2*Num_cell);
        for i = 1:Num_cell
            fprintf(fid,'%d %d\n',1,i-1);
        end
        
        fprintf(fid,'\n');
        fprintf(fid,'%s ','CELL_TYPES');
        fprintf(fid,'%d\n',Num_cell);
        for i = 1:Num_cell
            fprintf(fid,'%d\n',1);
        end
        
        
        fprintf(fid,'%s ','POINT_DATA');
        fprintf(fid,'%d\n',Num_cell);
        fprintf(fid,'%s\n','SCALARS wusRNALevel float');
        fprintf(fid,'%s\n','LOOKUP_TABLE default');
        for i = 1:Num_cell  
            fprintf(fid,'%f\n',wusR(i));
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'%s\n','SCALARS wusNucleiLevel float');
        fprintf(fid,'%s\n','LOOKUP_TABLE default');
        for i = 1:Num_cell  
            fprintf(fid,'%f\n',wn(i));
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'%s\n','SCALARS wusCytoplasmLevel float');
        fprintf(fid,'%s\n','LOOKUP_TABLE default');
        for i = 1:Num_cell  
            fprintf(fid,'%f\n',wc(i));
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'%s\n','SCALARS CLV3Level float');
        fprintf(fid,'%s\n','LOOKUP_TABLE default');
        for i = 1:Num_cell  
            fprintf(fid,'%f\n',clv(i));
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'%s\n','SCALARS ckLigandLevel float');
        fprintf(fid,'%s\n','LOOKUP_TABLE default');
        for i = 1:Num_cell  
            fprintf(fid,'%f\n',ck(i));
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'%s\n','SCALARS ckReceptorLevel float');
        fprintf(fid,'%s\n','LOOKUP_TABLE default');
        for i = 1:Num_cell  
            fprintf(fid,'%f\n',ckR(i));
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'%s\n','SCALARS ckSignalLevel float');
        fprintf(fid,'%s\n','LOOKUP_TABLE default');
        for i = 1:Num_cell  
            fprintf(fid,'%f\n',CK(i));
        end
        fprintf(fid,'\n');
        
        fclose(fid);
        
        
        fid = fopen(['wus_data_crosssection_T' num2str(it*dt) '.vtk'],'w+');
        fprintf(fid,'%s\n','# vtk DataFile Version 3.0');
        fprintf(fid,'%s\n','Points representing individual cells');
        fprintf(fid,'%s\n\n','ASCII');
        fprintf(fid,'%s\n','DATASET UNSTRUCTURED_GRID');
        fprintf(fid,'%s ','POINTS');
        fprintf(fid,'%d ',length(id_cs));
        fprintf(fid,'%s\n','float');
        for i = 1:length(id_cs)
            fprintf(fid,'%f %f %f\n',cell_xyz(id_cs(i),1),cell_xyz(id_cs(i),2),cell_xyz(id_cs(i),3));
        end
        
        fprintf(fid,'\n');
        fprintf(fid,'%s ','CELLS');
        fprintf(fid,'%d ',length(id_cs));
        fprintf(fid,'%d\n',2*length(id_cs));
        for i = 1:length(id_cs)
            fprintf(fid,'%d %d\n',1,i-1);
        end
        
        fprintf(fid,'\n');
        fprintf(fid,'%s ','CELL_TYPES');
        fprintf(fid,'%d\n',length(id_cs));
        for i = 1:length(id_cs)
            fprintf(fid,'%d\n',1);
        end
        
        
        fprintf(fid,'%s ','POINT_DATA');
        fprintf(fid,'%d\n',length(id_cs));
        fprintf(fid,'%s\n','SCALARS wusRNALevel float');
        fprintf(fid,'%s\n','LOOKUP_TABLE default');
        for i = 1:length(id_cs)  
            fprintf(fid,'%f\n',wusR(id_cs(i)));
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'%s\n','SCALARS wusNucleiLevel float');
        fprintf(fid,'%s\n','LOOKUP_TABLE default');
        for i = 1:length(id_cs)  
            fprintf(fid,'%f\n',wn(id_cs(i)));
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'%s\n','SCALARS wusCytoplasmLevel float');
        fprintf(fid,'%s\n','LOOKUP_TABLE default');
        for i = 1:length(id_cs)  
            fprintf(fid,'%f\n',wc(id_cs(i)));
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'%s\n','SCALARS CLV3Level float');
        fprintf(fid,'%s\n','LOOKUP_TABLE default');
        for i = 1:length(id_cs)  
            fprintf(fid,'%f\n',clv(id_cs(i)));
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'%s\n','SCALARS ckLigandLevel float');
        fprintf(fid,'%s\n','LOOKUP_TABLE default');
        for i = 1:length(id_cs)  
            fprintf(fid,'%f\n',ck(id_cs(i)));
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'%s\n','SCALARS ckReceptorLevel float');
        fprintf(fid,'%s\n','LOOKUP_TABLE default');
        for i = 1:length(id_cs) 
            fprintf(fid,'%f\n',ckR(id_cs(i)));
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'%s\n','SCALARS ckSignalLevel float');
        fprintf(fid,'%s\n','LOOKUP_TABLE default');
        for i = 1:length(id_cs)  
            fprintf(fid,'%f\n',CK(id_cs(i)));
        end
        fprintf(fid,'\n');
        
        clv_sc = clv_source(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3));
        clv_prod = clv_sc .* ( 2./( 1 + ((wthr1-wn)/klow).^(3*n)) ).*(wn<wthr1)+...
            clv_sc .* (2./( 1 + ((wn-wthr1)/(wthr2-wthr1)).^(3*n) )).*(wn>=wthr1);
        fprintf(fid,'%s\n','SCALARS CLV3_expression float');
        fprintf(fid,'%s\n','LOOKUP_TABLE default');
        for i = 1:length(id_cs)
            fprintf(fid,'%f\n',clv_prod(id_cs(i)));
        end
        fprintf(fid,'\n');
        
        wusR_prod = wusR_source(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3)) ...
            .* ( kcw1^(4*n)./(kcw1^(4*n)+clv.^(4*n)) );
        fprintf(fid,'%s\n','SCALARS WUS_expression float');
        fprintf(fid,'%s\n','LOOKUP_TABLE default');
        for i = 1:length(id_cs)
            fprintf(fid,'%f\n',wusR_prod(id_cs(i)));
        end
        fprintf(fid,'\n');
        
        fclose(fid);
    end
    
end





