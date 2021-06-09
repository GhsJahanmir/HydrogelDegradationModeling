close all 
clear all 
clc
time_exp=[     ]; %% the time point from experimental data
Swelling_mass_ratio_exp=[   ];    %%% the swelling ratio over dried mass in the corresponding time points
Swelling_ratio_exp=Swelling_mass_ratio_exp/Swelling_mass_ratio_exp(1);
%%%Initial Properties%%%%%%%%%
C_VS=35/100;%mg/ul 
C_SH=35/100;%mg/ul
DM_VS=0.09;    % degree of modification of VS polymer
DM_SH=0.09;    % degree of modification of SH polymer

GelMass_rel=30 ; %%%%mg or ul    relaxed state or volume of solution before gelation


V_gel=GelMass_rel;%ul  


Qm1=Swelling_mass_ratio_exp(1);


Mr_VS=162;% RU for one saccharide unit  
Mr_SH=162;  %
Mw_VS=40000;%Da   polymer MW
Mw_SH=40000;%Da   polymer MW
Vbar2_VS=.64;
Vbar2_SH=.64;


d_t=0.05;%day
end_t=13;%day


save('gel_expdata')

