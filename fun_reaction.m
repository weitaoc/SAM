function [fwusR,fwn,fwc,fclv,fck,fckR,fCK] = fun_reaction(cell_xyz,wusR,wn,wc,clv,ck,ckR,CK,wusR_source,clv_source,ck_source,ckR_source)
global EL WB

bio_parameters;

% WUS mRNA
wusR_prod = wusR_source(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3)) .* ( kcw1^(4*n)./(kcw1^(4*n)+clv.^(4*n)) ); % inhibition by clv
wusR_deg  = dw*wusR;
fwusR = wusR_prod - wusR_deg;

% WUS in nuclei
wn_prod   = rn*wusR;
% cytokinin improves the stability, CLV3 improves the stability,
% self-stabilization
wn_deg    = dwn.*( 1 + 2 ./ (1+(wn/kww).^(2*n)) ) .*wn;
rex       = 2*( ( rmin + ( rmax - rmin )/(1+(2*WB)^10) ) / (1+EL) );
%     rex       = 0.2+rex./(1+(r2d/20).^10); %differential export rate in different layers
wn_export = rex./( (1+(CK*WB/kckw1).^(2*n)) .* (1+(clv*EL/kcw2).^(5*n)) ) .* wn; % CLV3 reduces nuclear export
wn_import = rim .* wc;
fwn = wn_prod - wn_deg + wn_import - wn_export;

% WUS in cytoplasm
wc_prod = rc*wusR;
% cytokinin improves the stability, CLV3 stabilizes Wc,
% self-stabilization
wc_deg  = dwc./(1+(CK/kckw2).^(n)).*(0.05+0.95) .*wc;
fwc = wc_prod - wc_deg - wn_import + wn_export;

% CLV3
clv_sc = clv_source(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3));
clv_prod = clv_sc .* ( 2./( 1 + ((wthr1-wn)/klow).^(3*n)) ).*(wn<wthr1)+...
    clv_sc .* (2./( 1 + ((wn-wthr1)/(wthr2-wthr1)).^(3*n) )).*(wn>=wthr1); % activation by wus
clv_deg  = dc*clv;
fclv = clv_prod - clv_deg;

% cytokinin signaling
CK_prod = kon * ck .* ckR;
CK_deg  = dCK * CK;
fCK = CK_prod - CK_deg;

% cytokinin receptor
ckR_prod = ckR_source(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3));
ckR_deg  = dckR*ckR;
fckR = ckR_prod - ckR_deg - CK_prod;


% cytokinin ligand
ck_prod = ck_source(cell_xyz(:,1),cell_xyz(:,2),cell_xyz(:,3));
ck_deg  = dck*ck;
fck = ck_prod - ck_deg - CK_prod;
