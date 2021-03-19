% PRCC sensitivity analysis for LHS: inputs are rex/rex0, Dw/Dw0, dwn/dwn0,
% dwc/dwc0; outputs is (L3N,L2N,L1N)
clear all;

nsamp = 180; %sample size

cd sample_set1;
output_s1 = zeros(nsamp,3);
for ind_sample = 1:nsamp
    
    filename = ['s' num2str(ind_sample) '.mat'];
    load(filename);
    
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
    
    
    L3N = midline_wn(logical((midline_xyz(:,3)<7.5).*(midline_xyz(:,3)>5.5)));
    L2N = midline_wn(logical((midline_xyz(:,3)<8.5).*(midline_xyz(:,3)>7.5)));
    L1N = midline_wn((midline_xyz(:,3)>8.5));
    
    L3R = midline_wn(logical((midline_xyz(:,3)<7.5).*(midline_xyz(:,3)>5.5)))./midline_wc(logical((midline_xyz(:,3)<7.5).*(midline_xyz(:,3)>5.5)));
    L2R = midline_wn(logical((midline_xyz(:,3)<8.5).*(midline_xyz(:,3)>7.5)))./midline_wc(logical((midline_xyz(:,3)<8.5).*(midline_xyz(:,3)>7.5)));
    L1R = midline_wn((midline_xyz(:,3)>8.5))./midline_wc((midline_xyz(:,3)>8.5));
    
    output_s1(ind_sample,:) = [mean(L3N) mean(L2N) mean(L1N)];
end
cd ..


cd sample_set2;
output_s2 = zeros(nsamp,3);
for ind_sample = 1:nsamp
    
    filename = ['s' num2str(ind_sample) '.mat'];
    load(filename);
    
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
    
    
    L3N = midline_wn(logical((midline_xyz(:,3)<7.5).*(midline_xyz(:,3)>5.5)));
    L2N = midline_wn(logical((midline_xyz(:,3)<8.5).*(midline_xyz(:,3)>7.5)));
    L1N = midline_wn((midline_xyz(:,3)>8.5));
    
    L3R = midline_wn(logical((midline_xyz(:,3)<7.5).*(midline_xyz(:,3)>5.5)))./midline_wc(logical((midline_xyz(:,3)<7.5).*(midline_xyz(:,3)>5.5)));
    L2R = midline_wn(logical((midline_xyz(:,3)<8.5).*(midline_xyz(:,3)>7.5)))./midline_wc(logical((midline_xyz(:,3)<8.5).*(midline_xyz(:,3)>7.5)));
    L1R = midline_wn((midline_xyz(:,3)>8.5))./midline_wc((midline_xyz(:,3)>8.5));
    
    output_s2(ind_sample,:) = [mean(L3N) mean(L2N) mean(L1N)];
end
cd ..

% for ind_input = 1 : 5
%     
%     foldername = ['sample_setMixN' num2str(ind_input)];
%     cd(foldername);
% 
%     output_mix = zeros(nsamp,3);
%     for ind_sample = 1:nsamp
%         
%         filename = ['s' num2str(ind_sample) '.mat'];
%         load(filename);
%         
%         [ind] = find(cell_xyz(:,1).^2+cell_xyz(:,2).^2<=2.5);
%         midline_xyz = cell_xyz(ind,:);
%         midline_wusR = wusR(ind);
%         midline_wn = wn(ind);
%         midline_wc = wc(ind);
%         midline_clv = clv(ind);
%         midline_ck = ck(ind);
%         midline_ckR = ckR(ind);
%         midline_CK = CK(ind);
%         
%         [val,ind] = sort(midline_xyz(:,3),'ascend');
%         midline_xyz = midline_xyz(ind,:);
%         midline_wusR = midline_wusR(ind);
%         midline_wn = midline_wn(ind);
%         midline_wc = midline_wc(ind);
%         midline_clv = midline_clv(ind);
%         midline_ck = midline_ck(ind);
%         midline_ckR = midline_ckR(ind);
%         midline_CK = midline_CK(ind);
%         
%         
%         L3N = midline_wn(logical((midline_xyz(:,3)<7.5).*(midline_xyz(:,3)>5.5)));
%         L2N = midline_wn(logical((midline_xyz(:,3)<8.5).*(midline_xyz(:,3)>7.5)));
%         L1N = midline_wn((midline_xyz(:,3)>8.5));
%         
%         L3R = midline_wn(logical((midline_xyz(:,3)<7.5).*(midline_xyz(:,3)>5.5)))./midline_wc(logical((midline_xyz(:,3)<7.5).*(midline_xyz(:,3)>5.5)));
%         L2R = midline_wn(logical((midline_xyz(:,3)<8.5).*(midline_xyz(:,3)>7.5)))./midline_wc(logical((midline_xyz(:,3)<8.5).*(midline_xyz(:,3)>7.5)));
%         L1R = midline_wn((midline_xyz(:,3)>8.5))./midline_wc((midline_xyz(:,3)>8.5));
%         
%         output_mix(ind_sample,:) = [mean(L3N) mean(L2N) mean(L1N)];
%     end
%     varname = ['output_mix' num2str(ind_input)];
%     eval([varname '= output_mix;']);
%     cd ..
% end

for ind_input = 1 : 5
    
    foldername = ['sample_setMixNN' num2str(ind_input)];
    cd(foldername);

    output_mix = zeros(nsamp,3);
    for ind_sample = 1:nsamp
        
        filename = ['s' num2str(ind_sample) '.mat'];
        load(filename);
        
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
        
        
        L3N = midline_wn(logical((midline_xyz(:,3)<7.5).*(midline_xyz(:,3)>5.5)));
        L2N = midline_wn(logical((midline_xyz(:,3)<8.5).*(midline_xyz(:,3)>7.5)));
        L1N = midline_wn((midline_xyz(:,3)>8.5));
        
        L3R = midline_wn(logical((midline_xyz(:,3)<7.5).*(midline_xyz(:,3)>5.5)))./midline_wc(logical((midline_xyz(:,3)<7.5).*(midline_xyz(:,3)>5.5)));
        L2R = midline_wn(logical((midline_xyz(:,3)<8.5).*(midline_xyz(:,3)>7.5)))./midline_wc(logical((midline_xyz(:,3)<8.5).*(midline_xyz(:,3)>7.5)));
        L1R = midline_wn((midline_xyz(:,3)>8.5))./midline_wc((midline_xyz(:,3)>8.5));
        
        output_mix(ind_sample,:) = [mean(L3N) mean(L2N) mean(L1N)];
    end
    varname = ['output_mixx' num2str(ind_input)];
    eval([varname '= output_mix;']);
    cd ..
end


save('SS180.mat');

ind = find(isnan(output_s1(:,1)));
output_s1(ind,:) = [];
output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
output_mixx1(ind,:) = [];
output_mixx2(ind,:) = [];
output_mixx3(ind,:) = [];
output_mixx4(ind,:) = [];
output_mixx5(ind,:) = [];

ind = find(isnan(output_s2(:,1)));
output_s1(ind,:) = [];
output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
output_mixx1(ind,:) = [];
output_mixx2(ind,:) = [];
output_mixx3(ind,:) = [];
output_mixx4(ind,:) = [];
output_mixx5(ind,:) = [];

% ind = find(isnan(output_mix1(:,1)));
% output_s1(ind,:) = [];
% output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
% output_mixx1(ind,:) = [];
% output_mixx2(ind,:) = [];
% output_mixx3(ind,:) = [];
% output_mixx4(ind,:) = [];
% output_mixx5(ind,:) = [];
% 
% ind = find(isnan(output_mix2(:,1)));
% output_s1(ind,:) = [];
% output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
% output_mixx1(ind,:) = [];
% output_mixx2(ind,:) = [];
% output_mixx3(ind,:) = [];
% output_mixx4(ind,:) = [];
% output_mixx5(ind,:) = [];
% 
% ind = find(isnan(output_mix3(:,1)));
% output_s1(ind,:) = [];
% output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
% output_mixx1(ind,:) = [];
% output_mixx2(ind,:) = [];
% output_mixx3(ind,:) = [];
% output_mixx4(ind,:) = [];
% output_mixx5(ind,:) = [];
% 
% ind = find(isnan(output_mix4(:,1)));
% output_s1(ind,:) = [];
% output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
% output_mixx1(ind,:) = [];
% output_mixx2(ind,:) = [];
% output_mixx3(ind,:) = [];
% output_mixx4(ind,:) = [];
% output_mixx5(ind,:) = [];
% 
% ind = find(isnan(output_mix5(:,1)));
% output_s1(ind,:) = [];
% output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
% output_mixx1(ind,:) = [];
% output_mixx2(ind,:) = [];
% output_mixx3(ind,:) = [];
% output_mixx4(ind,:) = [];
% output_mixx5(ind,:) = [];

ind = find(isnan(output_mixx1(:,1)));
output_s1(ind,:) = [];
output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
output_mixx1(ind,:) = [];
output_mixx2(ind,:) = [];
output_mixx3(ind,:) = [];
output_mixx4(ind,:) = [];
output_mixx5(ind,:) = [];

ind = find(isnan(output_mixx2(:,1)));
output_s1(ind,:) = [];
output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
output_mixx1(ind,:) = [];
output_mixx2(ind,:) = [];
output_mixx3(ind,:) = [];
output_mixx4(ind,:) = [];
output_mixx5(ind,:) = [];

ind = find(isnan(output_mixx3(:,1)));
output_s1(ind,:) = [];
output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
output_mixx1(ind,:) = [];
output_mixx2(ind,:) = [];
output_mixx3(ind,:) = [];
output_mixx4(ind,:) = [];
output_mixx5(ind,:) = [];

ind = find(isnan(output_mixx4(:,1)));
output_s1(ind,:) = [];
output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
output_mixx1(ind,:) = [];
output_mixx2(ind,:) = [];
output_mixx3(ind,:) = [];
output_mixx4(ind,:) = [];
output_mixx5(ind,:) = [];

ind = find(isnan(output_mixx5(:,1)));
output_s1(ind,:) = [];
output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
output_mixx1(ind,:) = [];
output_mixx2(ind,:) = [];
output_mixx3(ind,:) = [];
output_mixx4(ind,:) = [];
output_mixx5(ind,:) = [];


level_threshold = 1000;
ind = find(output_s1(:,1)>level_threshold);
output_s1(ind,:) = [];
output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
output_mixx1(ind,:) = [];
output_mixx2(ind,:) = [];
output_mixx3(ind,:) = [];
output_mixx4(ind,:) = [];
output_mixx5(ind,:) = [];
ind = find(output_s2(:,1)>level_threshold);
output_s1(ind,:) = [];
output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
output_mixx1(ind,:) = [];
output_mixx2(ind,:) = [];
output_mixx3(ind,:) = [];
output_mixx4(ind,:) = [];
output_mixx5(ind,:) = [];
% ind = find(output_mix1(:,1)>level_threshold);
% output_s1(ind,:) = [];
% output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
% output_mixx1(ind,:) = [];
% output_mixx2(ind,:) = [];
% output_mixx3(ind,:) = [];
% output_mixx4(ind,:) = [];
% output_mixx5(ind,:) = [];
% ind = find(output_mix2(:,1)>level_threshold);
% output_s1(ind,:) = [];
% output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
% output_mixx1(ind,:) = [];
% output_mixx2(ind,:) = [];
% output_mixx3(ind,:) = [];
% output_mixx4(ind,:) = [];
% output_mixx5(ind,:) = [];
% ind = find(output_mix3(:,1)>level_threshold);
% output_s1(ind,:) = [];
% output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
% output_mixx1(ind,:) = [];
% output_mixx2(ind,:) = [];
% output_mixx3(ind,:) = [];
% output_mixx4(ind,:) = [];
% output_mixx5(ind,:) = [];
% ind = find(output_mix4(:,1)>level_threshold);
% output_s1(ind,:) = [];
% output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
% output_mixx1(ind,:) = [];
% output_mixx2(ind,:) = [];
% output_mixx3(ind,:) = [];
% output_mixx4(ind,:) = [];
% output_mixx5(ind,:) = [];
% ind = find(output_mix5(:,1)>level_threshold);
% output_s1(ind,:) = [];
% output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
% output_mixx1(ind,:) = [];
% output_mixx2(ind,:) = [];
% output_mixx3(ind,:) = [];
% output_mixx4(ind,:) = [];
% output_mixx5(ind,:) = [];
ind = find(output_mixx1(:,1)>level_threshold);
output_s1(ind,:) = [];
output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
output_mixx1(ind,:) = [];
output_mixx2(ind,:) = [];
output_mixx3(ind,:) = [];
output_mixx4(ind,:) = [];
output_mixx5(ind,:) = [];
ind = find(output_mixx2(:,1)>level_threshold);
output_s1(ind,:) = [];
output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
output_mixx1(ind,:) = [];
output_mixx2(ind,:) = [];
output_mixx3(ind,:) = [];
output_mixx4(ind,:) = [];
output_mixx5(ind,:) = [];
ind = find(output_mixx3(:,1)>level_threshold);
output_s1(ind,:) = [];
output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
output_mixx1(ind,:) = [];
output_mixx2(ind,:) = [];
output_mixx3(ind,:) = [];
output_mixx4(ind,:) = [];
output_mixx5(ind,:) = [];
ind = find(output_mixx4(:,1)>level_threshold);
output_s1(ind,:) = [];
output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
output_mixx1(ind,:) = [];
output_mixx2(ind,:) = [];
output_mixx3(ind,:) = [];
output_mixx4(ind,:) = [];
output_mixx5(ind,:) = [];
ind = find(output_mixx5(:,1)>level_threshold);
output_s1(ind,:) = [];
output_s2(ind,:) = [];
% output_mix1(ind,:) = [];
% output_mix2(ind,:) = [];
% output_mix3(ind,:) = [];
% output_mix4(ind,:) = [];
% output_mix5(ind,:) = [];
output_mixx1(ind,:) = [];
output_mixx2(ind,:) = [];
output_mixx3(ind,:) = [];
output_mixx4(ind,:) = [];
output_mixx5(ind,:) = [];

nsamp = size(output_s1,1);

% L3 nuclear level
for ind_input = 1:5
    
%     varname = ['output_mix' num2str(ind_input)];
%     eval(['output_mix=' varname ';']);
%     up = 1/(nsamp-1)*(output_s1(:,1)'*output_mix(:,1));
    Esquare = 1/nsamp*(output_s1(:,1)'*output_s2(:,1));
    V = 1/nsamp*(output_s1(:,1)'*output_s1(:,1))-(1/nsamp*sum(output_s1(:,1)))^2;
%     varname = ['L3SS' num2str(ind_input)];
%     eval([varname '=(up-Esquare)/V;']);
    
    varname = ['output_mixx' num2str(ind_input)];
    eval(['output_mix=' varname ';']);
    un = 1/(nsamp-1)*(output_s1(:,1)'*output_mix(:,1));
    varname = ['L3totSS' num2str(ind_input)];
    eval([varname '=1-(un-Esquare)/V;']);
    
end

% L2 nuclear level
for ind_input = 1:5
    
%     varname = ['output_mix' num2str(ind_input)];
%     eval(['output_mix=' varname ';']);
%     up = 1/(nsamp-1)*(output_s1(:,2)'*output_mix(:,2));
    Esquare = 1/nsamp*(output_s1(:,2)'*output_s2(:,2));
    V = 1/nsamp*(output_s1(:,2)'*output_s1(:,2))-(1/nsamp*sum(output_s1(:,2)))^2;
%     varname = ['L2SS' num2str(ind_input)];
%     eval([varname '=(up-Esquare)/V;']);
    
    varname = ['output_mixx' num2str(ind_input)];
    eval(['output_mix=' varname ';']);
    un = 1/(nsamp-1)*(output_s1(:,2)'*output_mix(:,2));
    varname = ['L2totSS' num2str(ind_input)];
    eval([varname '=1-(un-Esquare)/V;']);
end

% L1 nuclear level
for ind_input = 1:5
    
%     varname = ['output_mix' num2str(ind_input)];
%     eval(['output_mix=' varname ';']);
%     up = 1/(nsamp-1)*(output_s1(:,3)'*output_mix(:,3));
    Esquare = 1/nsamp*(output_s1(:,3)'*output_s2(:,3));
    V = 1/nsamp*(output_s1(:,3)'*output_s1(:,3))-(1/nsamp*sum(output_s1(:,3)))^2;
%     varname = ['L1SS' num2str(ind_input)];
%     eval([varname '=(up-Esquare)/V;']);
    
    varname = ['output_mixx' num2str(ind_input)];
    eval(['output_mix=' varname ';']);
    un = 1/(nsamp-1)*(output_s1(:,3)'*output_mix(:,3));
    varname = ['L1totSS' num2str(ind_input)];
    eval([varname '=1-(un-Esquare)/V;']);
    
end

% figure(3);
% bar([L3SS1 L3SS2 L3SS3 L3SS4 L3SS5])
% ylabel('nuclear WUS in L3');
% xticklabels({'synthesis','diffusion','nuclear degradation','cytoplasmic degradation','nuclear export'})
% set(gca,'FontSize',20,'LineWidth',2)
% ylim([-0.1 0.3])
% 
% 
% figure(2);
% bar([L2SS1 L2SS2 L2SS3 L2SS4 L2SS5])
% ylabel('nuclear WUS in L2');
% xticklabels({'synthesis','diffusion','nuclear degradation','cytoplasmic degradation','nuclear export'})
% set(gca,'FontSize',20,'LineWidth',2)
% ylim([-0.1 0.3])
% 
% 
% figure(1);
% bar([L1SS1 L1SS2 L1SS3 L1SS4 L1SS5])
% ylabel('nuclear WUS in L1');
% xticklabels({'synthesis','diffusion','nuclear degradation','cytoplasmic degradation','nuclear export'})
% set(gca,'FontSize',20,'LineWidth',2)
% ylim([-0.1 0.3])


figure(13);
bar([L3totSS1 L3totSS2 L3totSS3 L3totSS4 L3totSS5])
ylabel('nuclear WUS in L3');
xticklabels({'synthesis','diffusion','nuclear degradation','cytoplasmic degradation','nuclear export'})
set(gca,'FontSize',20,'LineWidth',2)
ylim([-1.5 1])


figure(12);
bar([L2totSS1 L2totSS2 L2totSS3 L2totSS4 L2totSS5])
ylabel('nuclear WUS in L2');
xticklabels({'synthesis','diffusion','nuclear degradation','cytoplasmic degradation','nuclear export'})
set(gca,'FontSize',20,'LineWidth',2)
ylim([-1.5 1])


figure(11);
bar([L1totSS1 L1totSS2 L1totSS3 L1totSS4 L1totSS5])
ylabel('nuclear WUS in L1');
xticklabels({'synthesis','diffusion','nuclear degradation','cytoplasmic degradation','nuclear export'})
set(gca,'FontSize',20,'LineWidth',2)
ylim([-1.5 1])
