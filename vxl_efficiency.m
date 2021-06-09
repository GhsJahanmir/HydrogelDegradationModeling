 
function [ve_real,eta]=vxl_efficiency(x_s,ve_ideal,V_gel,x_r)
V1=18;
X1_VS=0.46;
X1_SH=0.46;
X1=(X1_SH+X1_VS)/2;
%x_r=0.087225
%F_x=zeros(1,9);
%V=0.033;
%x0=ones(1,9)*xr;
%ve=93.15

ve_real=-(((log(1-x_s)+x_s+X1*x_s^2)/((x_s/x_r)^(1/3)-.5*(x_s/x_r)))/V1)*V_gel*10^-3*10^9;
eta=ve_real/ve_ideal
end