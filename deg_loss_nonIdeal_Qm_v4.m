function [time,Mole_current,Mass_current,P_0,X_s,ve_new,Qm_real,N_VSt,N_SHt,DM_VSt,DM_SHt,eta,t_c,Qm_ideal,Qv_real,ve_real1]=deg_loss_nonIdeal_Qm_v4(C_VS,C_SH,V_gel,DM_VS,DM_SH,Mw_VS,Mw_SH,Vbar2_VS,Vbar2_SH,k1,d_t,end_t,Qm1,Mr_VS,Mr_SH)
Vbar1=1.0069;
Vbar2=(Vbar2_VS+Vbar2_SH)/2; %ml/gr


Mass_VS=1000000*C_VS*V_gel/2;%ng
Mass_SH=1000000*C_SH*V_gel/2;%ng

N_C_VS=round(Mass_VS/Mw_VS); %nanomol
N_C_SH=round(Mass_SH/Mw_SH); %nanomol
N_VS=round((Mass_VS/Mr_VS*DM_VS)/N_C_VS);
N_SH=round((Mass_SH/Mr_SH*DM_SH)/N_C_SH);

Min_1=min(N_C_VS*N_VS,N_C_SH*N_SH);
Min=min((N_C_VS*N_C_SH),Min_1);


N_Ave_VS=round(Min/N_C_VS);
N_Ave_SH=round(Min/N_C_SH);



N_1=N_Ave_VS;
N_2=N_Ave_SH;

N_t=1+round(end_t/d_t);



if Min==N_C_SH*N_SH
    C_row=N_C_SH;
    C_column=N_C_VS;
    junction_node_r =N_2;% limit No for reacted points initially in one row
    junction_node_c =N_1;% limit No for reacted points initially in one column
    Nodes_0=round(rand(C_row,C_column));
    Nodes_i=Nodes_0;
end
if Min==N_C_VS*N_VS
    C_row=N_C_VS;
    C_column=N_C_SH;
    junction_node_r =N_1;% limit No for reacted points initially in one row
    junction_node_c =N_2;% limit No for reacted points initially in one column
    Nodes_0=round(rand(C_row,C_column));
    Nodes_i=Nodes_0;
end
if Min==N_C_VS*N_C_SH
    Min_3=min(N_1*N_C_VS,N_2*N_C_SH);
    if Min_3==N_1*N_C_VS
        C_row=N_C_VS;
        C_column=N_C_SH;
        junction_node_r =N_1;% limit No for reacted points initially
        junction_node_c =N_2;% limit No for reacted points initially
        Nodes_0=round(rand(C_row,C_column));
        Nodes_i=Nodes_0;
    else
        C_row=N_C_SH;
        C_column=N_C_VS;
        junction_node_r =N_2;% limit No for reacted points initially in one row
        junction_node_c =N_1;% limit No for reacted points initially in one column
        Nodes_0=round(rand(C_row,C_column));
        Nodes_i=Nodes_0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%obtaining initial network structure%%%%%%%%%%%%%%%%

In_r0=zeros(C_row,1);%intact bonds in row
In_c0=zeros(C_column,1);%intact bonds in column
In_r=zeros(C_row,1);%intact bonds in row
In_c=zeros(C_column,1);%intact bonds in column


for i=1:1:C_row
    
    In_r0(i)=sum(Nodes_0(i,:)==0);
    In_r(i)=In_r0(i);
end
for j=1:1:C_column
    In_c0(j)=sum(Nodes_0(:,j)==0);
    In_c(j)=In_c0(j);
end



for i=1:1:C_row
    for j=1:1:C_column
        if and(In_r(i)>junction_node_r,In_c(j)>junction_node_c)
            
            if Nodes_0(i,j)==0
                Nodes_i(i,j)=1;
                In_c(j)=In_c(j)-1;
                In_r(i)=In_r(i)-1;
                
            end
        end
        if and(In_c(j)<junction_node_c,In_r(i)<junction_node_r)
            
            if Nodes_0(i,j)==1
                Nodes_i(i,j)=0;
                In_c(j)=In_c(j)+1;
                In_r(i)=In_r(i)+1;
            end
        end
    end
end
for j=1:1:C_column
    for i=1:1:C_row
        if and(In_c(j)<junction_node_c,In_r(i)<junction_node_r)
            
            if Nodes_i(i,j)==1
                Nodes_i(i,j)=0;
                In_c(j)=In_c(j)+1;
                In_r(i)=In_r(i)+1;
            end
        end
        
        if and(In_c(j)>junction_node_c,Nodes_i(i,j)==0)
            
            Nodes_i(i,j)=1;
            In_c(j)=In_c(j)-1;
            In_r(i)=In_r(i)-1;
            
        end
        
    end
end
for i=1:1:C_row
    for j=1:1:C_column
        if and(In_r(i)>junction_node_r,Nodes_0(i,j)==0)
            
            Nodes_i(i,j)=1;
            In_c(j)=In_c(j)-1;
            In_r(i)=In_r(i)-1;
            
            
        end
        if and(In_c(j)<junction_node_c,In_r(i)<junction_node_r)
            
            if Nodes_0(i,j)==1
                Nodes_i(i,j)=0;
                In_c(j)=In_c(j)+1;
                In_r(i)=In_r(i)+1;
            end
        end
    end
end
for i=1:1:C_row
    for j=1:1:C_column
        
        
        if and(In_c(j)<junction_node_c,In_r(i)<junction_node_r)
            
            if Nodes_0(i,j)==1
                Nodes_i(i,j)=0;
                In_c(j)=In_c(j)+1;
                In_r(i)=In_r(i)+1;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  end   %%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%calculating efficiency of crosslinking and Ideal Swelling calculation %%%%%%%%%%
x_r=((C_VS+C_SH)/2)*Vbar2;
x_s1=Vbar2/(((Qm1-1)*Vbar1)+Vbar2);

for i=1:1:C_row
    if In_r(i)>0
        ve_r1(i)=In_r(i)-1;
    end
end

for j=1:1:C_column
    if In_c(j)>0
        ve_c1(j)=In_c(j)-1 ;
    end
end

ve_new_r1=sum(ve_r1(:));
ve_new_c1=sum(ve_c1(:));
ve_new1=ve_new_r1+ve_new_c1;
ve_ideal=ve_new1;
[ve_real1,eta]=vxl_efficiency(x_s1,ve_ideal,V_gel,x_r);


x_s_i=fsolve(@(x_s)Polymer_volFrac2(x_s,ve_ideal,V_gel,x_r),x_r);


Qm_ideal=(Vbar2*(1-x_s_i)+x_s_i*Vbar1)/(x_s_i*Vbar1);  %ideal Swelling ratio

%%%%%%%%%%%%%%%%%%%%%%end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%obtaining real network   %%%%%%%%%
if eta<1
    N_1_new=round((N_1-1)*eta+1); %real No. of crosslink points on VS polymer
    N_2_new=round((N_2-1)*eta+1); %real No. of crosslink points on SH polymer
else
    
    N_1_new=N_1; %real No. of crosslink points on VS polymer
    N_2_new=N_2; %real No. of crosslink points on SH polymer
end

%r=(N_1_new*N_C_VS)/(N_2_new*N_C_SH);
r=1; %% for hydrolysis
P=1/(r*(N_1_new-1)*(N_2_new-1))^.5;
t_c=-log(P)/k1; %disintegration time of gel


if Min==N_C_SH*N_SH
    C_row=N_C_SH;
    C_column=N_C_VS;
    junction_node_r =N_2_new;% limit No for reacted points initially in one row
    junction_node_c =N_1_new;% limit No for reacted points initially in one column
end
if Min==N_C_VS*N_VS
    C_row=N_C_VS;
    C_column=N_C_SH;
    junction_node_r =N_1_new;% limit No for reacted points initially in one row
    junction_node_c =N_2_new;% limit No for reacted points initially in one column
end
if Min==N_C_VS*N_C_SH
    Min_3=min(N_1*N_C_VS,N_2*N_C_SH);
    if Min_3==N_1*N_C_VS
        C_row=N_C_VS;
        C_column=N_C_SH;
        junction_node_r =N_1_new;% limit No for reacted points initially
        junction_node_c =N_2_new;% limit No for reacted points initially
    else
        C_row=N_C_SH;
        C_column=N_C_VS;
        junction_node_r =N_2_new;% limit No for reacted points initially in one row
        junction_node_c =N_1_new;% limit No for reacted points initially in one column
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%obtaining initial network structure%%%%%%%%%%%%%%%%

In_r_ideal=zeros(C_row,1);%intact bonds in row
In_c_ideal=zeros(C_column,1);%intact bonds in column
In_r_real=zeros(C_row,1);%intact bonds in row
In_c_real=zeros(C_column,1);%intact bonds in column
Nodes_r_i=ones(C_row,C_column); %real nodes of crosslink points
In_r_z=zeros(C_row,N_t);%intact bonds in row in different times
In_c_z=zeros(C_column,N_t);%intact bonds in column in different time

for i=1:1:C_row
    for j=1:1:C_column
        Nodes_r_i(i,j)=Nodes_i(i,j);
    end
end

for i=1:1:C_row
    In_r_ideal(i)=sum(Nodes_i(i,:)==0);
    In_r_real(i)=In_r_ideal(i);
end
for j=1:1:C_column
    In_c_ideal(j)=sum(Nodes_i(:,j)==0);
    In_c_real(j)=In_c_ideal(j);
end

for i=1:1:C_row
    for j=1:1:C_column
        if and(In_r_real(i)>junction_node_r,In_c_real(j)>junction_node_c)
            
            if Nodes_i(i,j)==0
                Nodes_r_i(i,j)=1;
                In_c_real(j)=In_c_real(j)-1;
                In_r_real(i)=In_r_real(i)-1;
                
            end
        end
        if and(In_c_real(j)<junction_node_c,In_r_real(i)<junction_node_r)
            
            if Nodes_i(i,j)==1
                Nodes_r_i(i,j)=0;
                In_c_real(j)=In_c_real(j)+1;
                In_r_real(i)=In_r_real(i)+1;
            end
        end
    end
end
for j=1:1:C_column
    for i=1:1:C_row
        if and(In_c_real(j)<junction_node_c,In_r_real(i)<junction_node_r)
            
            if Nodes_r_i(i,j)==1
                Nodes_r_i(i,j)=0;
                In_c_real(j)=In_c_real(j)+1;
                In_r_real(i)=In_r_real(i)+1;
            end
        end
        
        if and(In_c_real(j)>junction_node_c,Nodes_r_i(i,j)==0)
            
            Nodes_r_i(i,j)=1;
            In_c_real(j)=In_c_real(j)-1;
            In_r_real(i)=In_r_real(i)-1;
            
        end
        
    end
end
for i=1:1:C_row
    for j=1:1:C_column
        if and(In_r_real(i)>junction_node_r,Nodes_r_i(i,j)==0)
            
            Nodes_r_i(i,j)=1;
            In_c_real(j)=In_c_real(j)-1;
            In_r_real(i)=In_r_real(i)-1;
            
            
        end
        if and(In_c_real(j)<junction_node_c,In_r_real(i)<junction_node_r)
            
            if Nodes_r_i(i,j)==1
                Nodes_r_i(i,j)=0;
                In_c_real(j)=In_c_real(j)+1;
                In_r_real(i)=In_r_real(i)+1;
            end
        end
    end
end
for i=1:1:C_row
    for j=1:1:C_column
        
        
        if and(In_c_real(j)<junction_node_c,In_r_real(i)<junction_node_r)
            
            if Nodes_r_i(i,j)==1
                Nodes_r_i(i,j)=0;
                In_c_real(j)=In_c_real(j)+1;
                In_r_real(i)=In_r_real(i)+1;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  end  of obtaining real initial network %%%%%%%%%%%%%%%%


P_00=sum(sum(Nodes_i==0));%initial number of intact bonds in the network
P_real_00=sum(sum(Nodes_r_i==0));%initial number of intact bonds in the network
P_0=ones(1,N_t)*P_real_00;
P_frac=ones(1,N_t);
Mole_current_VS=ones(1,N_t)*N_C_VS;
Mole_current_SH=ones(1,N_t)*N_C_SH;
Mole_loss_VS=zeros(1,N_t);
Mole_loss_SH=zeros(1,N_t);
Mass_loss_VS=zeros(1,N_t);
Mass_loss_SH=zeros(1,N_t);

Mass_initial=ones(1,N_t)*(Mass_VS+Mass_SH);
Mass_current=Mass_initial;

Mole_initial=ones(1,N_t)*(N_C_VS+N_C_SH);
Mole_current=Mole_initial;

Mole_loss=zeros(1,N_t);
Mass_loss=zeros(1,N_t);

removed_row=zeros(1,N_t);
removed_column=zeros(1,N_t);

time=zeros(1,N_t);

ve=ones(1,N_t);
N_VSt=ones(1,N_t)*N_1_new;% No. of VS at time z on VS chain
N_SHt=ones(1,N_t)*N_2_new;% No. of SH at time z on SH chain
DM_VSt=ones(1,N_t)*N_1_new;% DM of VS at time z
DM_SHt=ones(1,N_t)*N_2_new;% DM of SH at time z

X_s=ones(1,N_t);
Qm_real=ones(1,N_t);
ve_new_r=zeros(1,N_t);
ve_new_c=zeros(1,N_t);
ve_new=zeros(1,N_t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %   life time of bonds%%%%%%%%%%%%%%%%%%%%%
Nodes_f=zeros(C_row,C_column);
t_nodes=zeros(C_row,C_column);
t_nodes_i=zeros(C_row,C_column);

for i=1:1:C_row
    for j=1:1:C_column
        %while( true )
        %eps=sigma.*randn(1,1)+mu;
        %eps=randn(1,1);
        % if (eps>0)&&(eps<1)
        %  if Nodes_f(i,j)==0
        %t_nodes_i(i,j)=-log(1-eps)/k1;
        %t_nodes(i,j)= t_nodes_i(i,j);
        % end
        % break
        %end
        % end
        
        Nodes_f(i,j)= Nodes_r_i(i,j);
        eps=0+rand(1,1)*(1-0);
        if Nodes_f(i,j)==0
            t_nodes_i(i,j)=-log(1-eps)/k1;
            t_nodes(i,j)= t_nodes_i(i,j);
        end
    end
end

%t_avg_0=mean(t_nodes(t_nodes~=0));
%t_avg=ones(1,N_t)*t_avg_0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  end of life time of bonds %%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%  start of iteration  %%%%%%%%%%%%%%%

for z=1:1:N_t
    
    t=(z-1)*d_t;
    for i=1:1:C_row
        for j=1:1:C_column
            if t_nodes(i,j)<=t
                if Nodes_f(i,j)==0
                    Nodes_f(i,j)=1;
                    In_r_real(i)=In_r_real(i)-1;
                    In_c_real(j)=In_c_real(j)-1;
                                                           
                    t_nodes(i,j)=0;
                    
                end
            end
        end
    end
    
    In_r_z(:,z)=In_r_real(:) ;
    In_c_z(:,z)=In_c_real(:);
    for i=1:1:C_row
        if In_r_real(i)>0
            ve_r(i)=In_r_real(i)-1;
        end
    end
    
    for j=1:1:C_column
        if In_c_real(j)>0
            ve_c(j)=In_c_real(j)-1 ;
        end
    end
    ve_new_r(z)=sum(ve_r(:));
    ve_new_c(z)=sum(ve_c(:));
    ve_new(z)=ve_new_r(z)+ve_new_c(z);
    for i=1:1:C_row
        if sum(Nodes_f(i,:))==C_column
            %or   if In_r_real(i)==0
            removed_row(z)=removed_row(z)+1;
        end
    end
    for j=1:1:C_column
        if sum(Nodes_f(:,j))==C_row
            %or   if In_c_real(j)==0
            removed_column(z)=removed_column(z)+1;
        end
    end
    
    P_0(z)=sum(sum(Nodes_f==0)); %number of intact bonds at each time points
    %  or P_0(z)=sum(In_r(:)) or P_0(z)=sum(In_c(:))
    P_frac(z)=P_0(z)/P_0(1);
    time(z)=t;
    
    
    
    if Min==N_C_VS*N_VS
        
        Mole_loss_VS(z)=(removed_row(z));
        Mole_loss_SH(z)=removed_column(z);
        Mole_current_VS(z)=(N_C_VS-Mole_loss_VS(z));
        Mole_current_SH(z)=(N_C_SH-Mole_loss_SH(z));
        Mole_loss(z)=Mole_loss_SH(z)+Mole_loss_VS(z);%nano mole
        Mole_current(z)=(Mole_initial(z))-Mole_loss(z);%nmol mole
        
        Mass_loss_VS(z)=(removed_row(z)*Mw_VS)*10^-6;
        Mass_loss_SH(z)=(removed_column(z)*Mw_SH)*10^-6;
        Mass_loss(z)=Mass_loss_SH(z)+Mass_loss_VS(z); %ng
        Mass_current(z)=(Mass_initial(z)*10^-6)-Mass_loss(z);
        
        ve(1)=4*(N_C_VS-1)*(N_1_new-1)-(N_C_VS-2)*(N_1_new-1)-(N_C_VS-1)*(N_1_new-2);%number of initial subchains
        N_VSt(z)=P_0(z)/Mole_current_VS(z);
        N_SHt(z)=P_0(z)/Mole_current_SH(z);
        DM_VSt(z)=100*N_VSt(z)*Mr_VS/Mw_VS;  %    %DM  of VS at time z
        DM_SHt(z)=100*N_SHt(z)*Mr_SH/Mw_SH;  %    %DM of VS at time z
        
        
        ve_1=4*(Mole_current_VS(z)-1)*(N_VSt(z)-1)-(Mole_current_VS(z)-2)*(N_VSt(z)-1)-(Mole_current_VS(z)-1)*(N_VSt(z)-2);
        ve_2=4*(Mole_current_SH(z)-1)*(N_SHt(z)-1)-(Mole_current_SH(z)-2)*(N_SHt(z)-1)-(Mole_current_SH(z)-1)*(N_SHt(z)-2);
        ve(z)=(ve_1+ve_2)/2; %number of subchains at time z
    end
    
    if Min==N_C_SH*N_SH
        Mole_loss_SH(z)=(removed_row(z));
        Mole_loss_VS(z)=removed_column(z);
        Mole_current_SH(z)=(N_C_SH-Mole_loss_SH(z));
        Mole_current_VS(z)=(N_C_VS-Mole_loss_VS(z));
        Mole_loss(z)=Mole_loss_SH(z)+Mole_loss_VS(z);%nano mole
        Mole_current(z)=(Mole_initial(z))-Mole_loss(z);%nmol mole
        
        Mass_loss_SH(z)=(removed_row(z)*Mw_SH)*10^-6;
        Mass_loss_VS(z)=(removed_column(z)*Mw_VS)*10^-6;
        Mass_loss(z)=Mass_loss_SH(z)+Mass_loss_VS(z); %ng
        Mass_current(z)=(Mass_initial(z)*10^-6)-Mass_loss(z);
        
        ve(1)=4*(N_C_SH-1)*(N_2_new-1)-(N_C_SH-2)*(N_2_new-1)-(N_C_SH-1)*(N_2_new-2);
        
        
        N_VSt(z)=P_0(z)/Mole_current_VS(z);
        N_SHt(z)=P_0(z)/Mole_current_SH(z);
        DM_VSt(z)=100*N_VSt(z)*Mr_VS/Mw_VS;  %    %DM  of VS at time z
        DM_SHt(z)=100*N_SHt(z)*Mr_SH/Mw_SH;  %    %DM of VS at time z
        
        ve_1=4*(Mole_current_SH(z)-1)*(N_SHt(z)-1)-(Mole_current_SH(z)-2)*(N_SHt(z)-1)-(Mole_current_SH(z)-1)*(N_SHt(z)-2);
        ve_2=4*(Mole_current_VS(z)-1)*(N_VSt(z)-1)-(Mole_current_VS(z)-2)*(N_VSt(z)-1)-(Mole_current_VS(z)-1)*(N_VSt(z)-2);
        ve(z)=(ve_1+ve_2)/2;  %number of subchains at time z
        
    end
    if Min==N_C_VS*N_C_SH
        
        Min_3=min(N_1*N_C_VS,N_2*N_C_SH);
        
        if Min_3==N_1*N_C_VS
            Mole_loss_VS(z)=(removed_row(z));
            Mole_loss_SH(z)=removed_column(z);
            Mole_current_VS(z)=(N_C_VS-Mole_loss_VS(z));
            Mole_current_SH(z)=(N_C_SH-Mole_loss_SH(z));
            Mole_loss(z)=Mole_loss_SH(z)+Mole_loss_VS(z);%nano mole
            Mole_current(z)=(Mole_initial(z))-Mole_loss(z);%nmol mole
            
            Mass_loss_VS(z)=(removed_row(z)*Mw_VS)*10^-6;
            Mass_loss_SH(z)=(removed_column(z)*Mw_SH)*10^-6;
            Mass_loss(z)=Mass_loss_SH(z)+Mass_loss_VS(z); %ng
            Mass_current(z)=(Mass_initial(z)*10^-6)-Mass_loss(z);
            
            ve(1)=4*(N_C_VS-1)*(N_1_new-1)-(N_C_VS-2)*(N_1_new-1)-(N_C_VS-1)*(N_1_new-2);%number of initial subchains
            N_VSt(z)=P_0(z)/Mole_current_VS(z);
            N_SHt(z)=P_0(z)/Mole_current_SH(z);
            DM_VSt(z)=100*N_VSt(z)*Mr_VS/Mw_VS;  %    %DM  of VS at time z
            DM_SHt(z)=100*N_SHt(z)*Mr_SH/Mw_SH;  %    %DM of VS at time z
            
            ve_1=4*(Mole_current_VS(z)-1)*(N_VSt(z)-1)-(Mole_current_VS(z)-2)*(N_VSt(z)-1)-(Mole_current_VS(z)-1)*(N_VSt(z)-2);
            ve_2=4*(Mole_current_SH(z)-1)*(N_SHt(z)-1)-(Mole_current_SH(z)-2)*(N_SHt(z)-1)-(Mole_current_SH(z)-1)*(N_SHt(z)-2);
            ve(z)=(ve_1+ve_2)/2; %number of subchains at time z
        else
            
            
            Mole_loss_SH(z)=(removed_row(z));
            Mole_loss_VS(z)=removed_column(z);
            Mole_current_SH(z)=(N_C_SH-Mole_loss_SH(z));
            Mole_current_VS(z)=(N_C_VS-Mole_loss_VS(z));
            Mole_loss(z)=Mole_loss_SH(z)+Mole_loss_VS(z);%nano mole
            Mole_current(z)=(Mole_initial(z))-Mole_loss(z);%nmol mole
            
            Mass_loss_SH(z)=(removed_row(z)*Mw_SH)*10^-6;
            Mass_loss_VS(z)=(removed_column(z)*Mw_VS)*10^-6;
            Mass_loss(z)=Mass_loss_SH(z)+Mass_loss_VS(z); %ng
            Mass_current(z)=(Mass_initial(z)*10^-6)-Mass_loss(z);
            
            ve(1)=4*(N_C_SH-1)*(N_2_new-1)-(N_C_SH-2)*(N_2_new-1)-(N_C_SH-1)*(N_2_new-2);
            
            
            N_VSt(z)=P_0(z)/Mole_current_VS(z);
            N_SHt(z)=P_0(z)/Mole_current_SH(z);
            DM_VSt(z)=100*N_VSt(z)*Mr_VS/Mw_VS;  %    %DM  of VS at time z
            DM_SHt(z)=100*N_SHt(z)*Mr_SH/Mw_SH;  %    %DM of VS at time z
            
            ve_1=4*(Mole_current_SH(z)-1)*(N_SHt(z)-1)-(Mole_current_SH(z)-2)*(N_SHt(z)-1)-(Mole_current_SH(z)-1)*(N_SHt(z)-2);
            ve_2=4*(Mole_current_VS(z)-1)*(N_VSt(z)-1)-(Mole_current_VS(z)-2)*(N_VSt(z)-1)-(Mole_current_VS(z)-1)*(N_VSt(z)-2);
            ve(z)=(ve_1+ve_2)/2;  %number of subchains at time z
            
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Swelling calculation   %%%%%%%%%%%%%%%%
    
    x_s=fsolve(@(x_s)Polymer_volFrac2(x_s,ve_new(z),V_gel,x_r),x_r);
    
    X_s(z)=x_s  ;
    q=(Vbar2*(1-x_s)+x_s*Vbar1)/(x_s*Vbar1);
    Qm_real(z)=q ; %Swelling ratio
    Qv_real(z)=x_r/x_s;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    
    
    
    if and(N_SHt(z)<1,N_VSt(z)<1)
        break
    end
    
end



end




