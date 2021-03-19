% plot along the midline x=0, y=0
clear all
load data_T270.mat;
[ind] = find(cell_xyz(:,1).^2+cell_xyz(:,2).^2<=2.5);
midline_xyz = cell_xyz(ind,:);
midline_wusR = wusR(ind);
midline_wn = wn(ind);
midline_wc = wc(ind);
midline_clv = clv(ind);
midline_ck = ck(ind);
midline_ckR = ckR(ind);
midline_CK = CK(ind);

[val,ind] = sort(midline_xyz(:,3),'ascend');
midline_xyz = midline_xyz(ind,:);
midline_wusR = midline_wusR(ind);
midline_wn = midline_wn(ind);
midline_wc = midline_wc(ind);
midline_clv = midline_clv(ind);
midline_ck = midline_ck(ind);
midline_ckR = midline_ckR(ind);
midline_CK = midline_CK(ind);

figure;
subplot(2,2,1)
plot(midline_xyz(:,3),midline_wusR,'o-','LineWidth',2);
subplot(2,2,2)
plot(midline_xyz(:,3),midline_clv,'o-','LineWidth',2);
subplot(2,2,3); 
yyaxis left
plot(midline_xyz(:,3),midline_wn,'o-')
yyaxis right
plot(midline_xyz(:,3),midline_wc,'*-');
subplot(2,2,4)
plot(midline_xyz(:,3),midline_CK,'o-','LineWidth',2);


figure;%expression domain
clv_sc = clv_source(midline_xyz(:,1),midline_xyz(:,2),midline_xyz(:,3));
clv_prod = clv_sc .* ( 2./( 1 + ((wthr1-midline_wn)/klow).^(3*n)) ).*(midline_wn<wthr1)+...
    clv_sc .* (2./( 1 + ((midline_wn-wthr1)/(wthr2-wthr1)).^(3*n) )).*(midline_wn>=wthr1);
subplot(1,2,1);
plot(midline_xyz(:,3),clv_prod,'*-');
wusR_prod = wusR_source(midline_xyz(:,1),midline_xyz(:,2),midline_xyz(:,3)) .* ( kcw1^(4*n)./(kcw1^(4*n)+midline_clv.^(4*n)) );
subplot(1,2,2);
plot(midline_xyz(:,3),wusR_prod,'*-');


L3N = midline_wn(logical((midline_xyz(:,3)<7.5).*(midline_xyz(:,3)>6.5)));
L2N = midline_wn(logical((midline_xyz(:,3)<8.5).*(midline_xyz(:,3)>7.5)));
L1N = midline_wn((midline_xyz(:,3)>8.5));
figure(10);hold on;
errorbar([mean(L3N) mean(L2N) mean(L1N)],[std(L3N) std(L2N) std(L1N)]);
disp(['Wn L3/L1 is ' num2str(mean(L3N)/mean(L1N))]);
disp(['Wn L2/L1 is ' num2str(mean(L2N)/mean(L1N))]);


L3C = midline_wc(logical((midline_xyz(:,3)<7.5).*(midline_xyz(:,3)>6.5)));
L2C = midline_wc(logical((midline_xyz(:,3)<8.5).*(midline_xyz(:,3)>7.5)));
L1C = midline_wc((midline_xyz(:,3)>8.5));
L3R = midline_wn(logical((midline_xyz(:,3)<7.5).*(midline_xyz(:,3)>6.5)))./midline_wc(logical((midline_xyz(:,3)<7.5).*(midline_xyz(:,3)>6.5)));
L2R = midline_wn(logical((midline_xyz(:,3)<8.5).*(midline_xyz(:,3)>7.5)))./midline_wc(logical((midline_xyz(:,3)<8.5).*(midline_xyz(:,3)>7.5)));
L1R = midline_wn((midline_xyz(:,3)>8.5))./midline_wc((midline_xyz(:,3)>8.5));
disp(['L3 NC ratio is ' num2str(mean(L3R))]);
disp(['L2 NC ratio is ' num2str(mean(L2R))]);
disp(['L1 NC ratio is ' num2str(mean(L1R))]);

figure(12);%show CLV3 regulations
trans_reg = kcw1^(4*n)./(kcw1^(4*n)+midline_clv.^(4*n));
block_exp = rmin ./ (1+(midline_clv/kcw2).^(5*n));
stab_cyt = 1./(1+(midline_clv/kcw3).^(5*n));
subplot(1,3,1)
plot(midline_xyz(:,3),trans_reg,'*-');
subplot(1,3,2)
plot(midline_xyz(:,3),block_exp,'*-');
subplot(1,3,3)
plot(midline_xyz(:,3),stab_cyt,'*-');


save('wt_L3N.txt','L3N','-ascii');
save('wt_L2N.txt','L2N','-ascii');
save('wt_L1N.txt','L1N','-ascii');
save('wt_L3C.txt','L3C','-ascii');
save('wt_L2C.txt','L2C','-ascii');
save('wt_L1C.txt','L2C','-ascii');
