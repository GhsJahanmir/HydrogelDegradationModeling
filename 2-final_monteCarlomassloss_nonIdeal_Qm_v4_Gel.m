close all
clear all
clc

load gel_expdata

N_t=1+round(end_t/d_t);
filename='Matlab_results_Gel.xlsx';
ERow=1;
n=100;
nn=6;

Mole_current_n=ones(n,N_t);    Mass_current_n=ones(n,N_t);     ve_new_n=ones(n,N_t);     X_s_n=ones(n,N_t);    P_0_n=ones(n,N_t);    Qm_real_n=ones(n,N_t);
N_VSt_n=ones(n,N_t); N_SHt_n=ones(n,N_t);  t_c_n=ones(n,1);     DM_VSt_n=ones(n,N_t); DM_SHt_n=ones(n,N_t);
eta_n=ones(n,1); Qm_ideal_n=ones(n,1); Qv_real_n=ones(n,N_t); ve_real1_n=ones(n,N_t);

Mole_current_ave=ones(nn,N_t); Mass_current_ave=ones(nn,N_t);  ve_new_ave=ones(nn,N_t);  X_s_ave=ones(nn,N_t); P_0_ave=ones(nn,N_t); Qm_real_ave=ones(nn,N_t);
N_VSt_ave=ones(nn,N_t); N_SHt_ave=ones(nn,N_t); t_c_ave=ones(nn,1);   DM_VSt_ave=ones(nn,N_t); DM_SHt_ave=ones(nn,N_t);
eta_ave=ones(nn,1);  Qm_ideal_ave=ones(nn,1); Qv_real_ave=ones(nn,N_t); ve_real1_ave=ones(nn,N_t);
Qm_real_ave_ratio=zeros(nn,N_t);
for j=1:1:nn
    
    k=0.1;   % day^-1  hydrolysis rate constant   %%%optimizing parameter
    k1=k+0.0025*(j-1);
    
    for i=1:1:n
        [time,Mole_current,Mass_current,P_0,X_s,ve_new,Qm_real,N_VSt,N_SHt,DM_VSt,DM_SHt,eta,t_c,Qm_ideal,Qv_real,ve_real1]=deg_loss_nonIdeal_Qm_v4(C_VS,C_SH,GelMass_rel,DM_VS,DM_SH,Mw_VS,Mw_SH,Vbar2_VS,Vbar2_SH,k1,d_t,end_t,Qm1,Mr_VS,Mr_SH); % Whatever function you call...
        
        Mole_current_n(i,:)=Mole_current;
        Mass_current_n(i,:)=Mass_current;
        ve_new_n(i,:)=ve_new;
        X_s_n(i,:)=X_s;
        P_0_n(i,:)=P_0;
        Qm_real_n(i,:)=Qm_real;
        N_VSt_n(i,:)=N_VSt;
        N_SHt_n(i,:)=N_SHt;
        DM_VSt_n(i,:)=DM_VSt;
        DM_SHt_n(i,:)=DM_SHt;
        eta_n(i,1)=eta;
        t_c_n(i,:)=t_c;
        Qm_ideal_n(i,:)=Qm_ideal;
        Qv_real_n(i,:)=Qv_real;
        ve_real1_n(i,:)=ve_real1;
        
        
        
    end
    
    
    cellReference_t=sprintf('A%d',j+2);
    
    xlswrite(filename,k1,1,cellReference_t);
    xlswrite(filename,k1,2,cellReference_t);
    xlswrite(filename,k1,3,cellReference_t);
    xlswrite(filename,k1,4,cellReference_t);
    xlswrite(filename,k1,5,cellReference_t);
    xlswrite(filename,k1,6,cellReference_t);
    xlswrite(filename,k1,7,cellReference_t);
    xlswrite(filename,k1,8,cellReference_t);
    xlswrite(filename,k1,9,cellReference_t);
    xlswrite(filename,k1,10,cellReference_t);
    xlswrite(filename,k1,11,cellReference_t);
    xlswrite(filename,k1,12,cellReference_t);
    xlswrite(filename,k1,13,cellReference_t);
    
    Mole_current_ave(j,:)=mean(Mole_current_n,1);
    Mass_current_ave(j,:)=mean(Mass_current_n,1);
    ve_new_ave(j,:)=mean(ve_new_n,1);
    X_s_ave(j,:)=mean(X_s_n,1);
    P_0_ave(j,:)=mean(P_0_n,1);
    Qm_real_ave(j,:)=mean(Qm_real_n,1);
    N_SHt_ave(j,:)=mean(N_SHt_n,1);
    N_VSt_ave(j,:)=mean(N_VSt_n,1);
    DM_SHt_ave(j,:)=mean(DM_SHt_n,1);
    DM_VSt_ave(j,:)=mean(DM_VSt_n,1);
    eta_ave(j,1)=mean(eta_n,1);
    t_c_ave(j,:)=mean(t_c_n,1);
    Qm_ideal_ave(j,:)=mean(Qm_ideal_n,1);
    Qv_real_ave(j,:)=mean(Qv_real_n,1);
    ve_real1_ave(j,:)=mean(ve_real1_n,1);
    
    
    if eta_ave(j,1)>1
        %when there is a deviation between exp and model initial Q, initial expQ< ideal Q of model
        Qm_real_ave_ratio(j,:)=Qm_real_ave(j,:)/Qm_real_ave(j,1);
        
        swelling_ratio_eq=pchip((time),Qm_real_ave_ratio(j,:),time_exp);  %%% find the exact fraction in corresponding experimental time points
        
        differ=swelling_ratio_eq-Swelling_ratio_exp;
        %%%%%%%  finding error between exp data and simulations%%%%%%%%%%%%%%%%%%%
        %%%%%%%  finding error between exp data and simulations%%%%%%%%%%%%%%%%%%%
        Rsquared(j)=1-(sum((swelling_ratio_eq-Swelling_ratio_exp).^2)/sum((Swelling_ratio_exp-mean(Swelling_ratio_exp)).^2))
        
    else
        swelling_ratio_eq(j,:)=pchip((time),Qm_real_ave(j,:),time_exp);
        differ(j,:)=swelling_ratio_eq(j,:)-Swelling_mass_ratio_exp;
        Sq_error1(j)=sum(differ(j,:).^2);
        Rsquared1(j)=1-(sum((swelling_ratio_eq(j,:)-Swelling_mass_ratio_exp).^2)/sum((Swelling_mass_ratio_exp-mean(Swelling_mass_ratio_exp)).^2))
    end
    %%%%%%% finding error between exp data and simulations%%%%%%%%%%%%%%%%%%%
    %%%%%%%  finding error between exp data and simulations%%%%%%%%%%%%%%%%%%%
    %Rsquared=1-(sum((swelling_ratio_eq-Swelling_ratio_exp).^2)/sum((Swelling_ratio_exp-mean(Swelling_ratio_exp)).^2))
    
    
    
end
cellReference_t=sprintf('B%d',1);

xlswrite(filename,time,1,cellReference_t);%1 is first excel sheet name
xlswrite(filename,time,2,cellReference_t);%2 is excel sheet name
xlswrite(filename,time,3,cellReference_t);%3 is excel sheet name
xlswrite(filename,time,4,cellReference_t);%4 is excel sheet name
xlswrite(filename,time,5,cellReference_t);%6 is excel sheet name
xlswrite(filename,time,6,cellReference_t);%7 is excel sheet name
xlswrite(filename,time,7,cellReference_t);%7 is excel sheet name
xlswrite(filename,time,8,cellReference_t);%7 is excel sheet name
xlswrite(filename,time,9,cellReference_t);%7 is excel sheet name
xlswrite(filename,time,10,cellReference_t);%7 is excel sheet name
xlswrite(filename,time,11,cellReference_t);%7 is excel sheet name
xlswrite(filename,time,12,cellReference_t);%7 is excel sheet name


index=sprintf('B%d',ERow+2);

xlswrite(filename,Mole_current_ave,1,index);%1  is excel sheet name put mole in
xlswrite(filename,Mass_current_ave,2,index);%2 is excel sheet name put Mass in
xlswrite(filename,P_0_ave,3,index);%3 is excel sheet name put ve(subchains) in
xlswrite(filename,X_s_ave,4,index);%4 is excel sheet name put ve(subchains) in
xlswrite(filename,ve_new_ave,5,index);%6 is excel sheet name put ve_new in
xlswrite(filename,Qm_real_ave,6,index);
xlswrite(filename,Qv_real_ave,7,index);
xlswrite(filename,N_VSt_ave,8,index);
xlswrite(filename,N_SHt_ave,9,index);
xlswrite(filename,DM_VSt_ave,10,index);
xlswrite(filename,DM_SHt_ave,11,index);
xlswrite(filename,Qm_ideal_ave,12,index);%7 is excel sheet name
xlswrite(filename,eta_ave,13,index);%5 is excel sheet name put ve_new in
xlswrite(filename,t_c_ave,14,index);
xlswrite(filename,ve_real1,15,index);



figure(1)
plot(time_exp,swelling_ratio_eq,'g-')
hold on

plot(time_exp,Swelling_mass_ratio_exp','ro');
hold off
title(['R^2 = ' num2str(Rsquared1)])
ylabel('Released Fraction','FontSize',12);
xlabel('time(day)','FontSize',12);
legend({'Model result','Experimental data'},'FontSize',12,'TextColor','black','Location','southeast');


save('Gel_FinalData')

winopen(filename)