clear all
global SITENAMES
SITENAMES={'950','970','997','1007','1060'}
A_Mon_coop_vec10=[1.1];
A_Dim_coop_vec10=[1.1];



B_Mon_coop_vec10=[50];
B_Dim_coop_vec10=[50];

Kp_Time10 = [0.2];

TimeMaxBind_vec10 = [10];
trans_time_vec10 = [4];

konc=0.1;

WUS=(0.008:0.01:0.4080);

for pm_i10 = 1:length(TimeMaxBind_vec10)
    for adc10 = 1:length(A_Dim_coop_vec10)
        for amc10 = 1:length(A_Mon_coop_vec10)
            for bdc10 = 1:length(B_Dim_coop_vec10)
                for bmc10 = 1:length(B_Mon_coop_vec10)
                    for ti10 = 1:length(trans_time_vec10)
                        for kpi10 = 1:length(Kp_Time10)
                            
                            
                            a_mon_coop10 = A_Mon_coop_vec10(amc10);
                            a_dim_coop10 = A_Dim_coop_vec10(adc10);
                            b_mon_coop10 = B_Mon_coop_vec10(bmc10);
                            b_dim_coop10 = B_Dim_coop_vec10(bdc10);
                            
                            kp10 = Kp_Time10(kpi10);
                            TimeMaxBind10 = TimeMaxBind_vec10(pm_i10);
                            trans_time10 =trans_time_vec10(ti10);
                            
                            
                            
                            figure
                            
                            str = ['WT_pm',num2str(TimeMaxBind10),'kp',num2str(kp10),'tr_t',num2str(trans_time10),'amc',num2str(a_mon_coop10),'adc',num2str(a_dim_coop10),'bmc',num2str(b_mon_coop10),'bdc',num2str(b_dim_coop10),'.mat'];
                            load(str)
                            
                            WTmRNA = (1/nr_sim)*summRNA;
                            
                            plot(WUS,WTmRNA,'Color',[0 0 0 ],'LineWidth',3)
                            
                            
                            hold on
                            
                            str = ['970M_pm',num2str(TimeMaxBind10),'kp',num2str(kp10),'tr_t',num2str(trans_time10),'amc',num2str(a_mon_coop10),'adc',num2str(a_dim_coop10),'bmc',num2str(b_mon_coop10),'bdc',num2str(b_dim_coop10),'.mat'];
                            
                            load(str)
                            M970mRNA = (1/nr_sim)*summRNA;
                            
                            plot(WUS,M970mRNA,'Color',[0 0 1 ],'LineWidth',3)
                            hold on
                            
                            str = ['950M_pm',num2str(TimeMaxBind10),'kp',num2str(kp10),'tr_t',num2str(trans_time10),'amc',num2str(a_mon_coop10),'adc',num2str(a_dim_coop10),'bmc',num2str(b_mon_coop10),'bdc',num2str(b_dim_coop10),'.mat'];
                            load(str)
                            M950mRNA = (1/nr_sim)*summRNA;
                            
                            plot(WUS,M950mRNA,'Color',[1 0 1 ],'LineWidth',3)
                            hold on
                            
                            
                            %triple mutant
                            
                            str = ['DM_pm',num2str(TimeMaxBind10),'kp',num2str(kp10),'tr_t',num2str(trans_time10),'amc',num2str(a_mon_coop10),'adc',num2str(a_dim_coop10),'bmc',num2str(b_mon_coop10),'bdc',num2str(b_dim_coop10),'.mat'];
                            load(str)
                            triplemRNA = (1/nr_sim)*summRNA;
                            
                            plot(WUS,triplemRNA,'Color',[0.5 0 0.6],'LineWidth',3)
                            hold on
                            
                            
                            legend({'WT','970M','950M','DM'},'FontSize',10)
                            
                            ylabel('mRNA')
                            xlabel('WUS')
                            set(gca,'FontSize',20)
                            xlim([0 0.2])
                            
                            
                            currentFolder = pwd;
                            
                            fname_fig = [currentFolder, '\figures\fig_files\'];
                            str1= ['multiple',num2str(TimeMaxBind10),'kp',num2str(kp10),'mc',num2str(b_mon_coop10),'dc',num2str(b_dim_coop10)];
                            strfig = [str1,'.fig'];
                            saveas(gca,strcat(fname_fig, strfig));
                            
                            fname_png = [currentFolder, '\figures\png_files\'];
                            strpng = [str1,'.png'];
                            saveas(gca,strcat(fname_png, strpng));
                            
                            
                            fname_eps =  [currentFolder, '\figures\eps_files\'];
                            streps = [str1,'.eps'];
                            saveas(gca,strcat(fname_eps, streps));
                            
                            fname_pdf =  [currentFolder, '\figures\pdf_files\'];
                            strpdf = [str1,'.pdf'];
                            saveas(gca,strcat(fname_pdf, strpdf));
                            
                        end
                    end
                end
            end
        end
    end
end