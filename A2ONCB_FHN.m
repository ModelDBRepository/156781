function A2ONCB_FHN
%%% In order to understand the slow platau-like oscillation arising in rd
%%% AII after LP application, we hypothesized that there is a slow Ca2+
%%% mediated oscillation in the coupled ONCB. We model the Ca oscillation
%%% as a relaxation oscillation based on FitzHugh-Nagumo model.


clc
clear all

global epsilon A vz0 vz1 vz2 vz3 vz4 w0 gamma v0 exp_amp exp_thresh exp_steep tau_w
global delta_tau v_tau v_xi delta_tauM v_tauM v_xiM

global gjONCB ie ie_tonic ie_pulse t_min t_max

global numAII armsection R capacitance            
global gbar_na gbar_M gbar_A G_gap Gpas epas ena ek      
global vhalfh_na vhalfm_na kh_na km_na vhalfm_M km_M vhalfh_A vhalfm_A vhalfw_A kh_A km_A kw_A f
global mtau_na htau_na mtau_M mtau_A 
global tau capacitanceONCB 
global t_switch t_width gbar_M_IS_1 gbar_M_IS_2 SA_hand


%%%%%%%% -------------------- Switches ---------------- %%%%%%%%%

numAII = 1;
ONCB = 10;%10;
ie = 0;
% current injection during the duration t_min to t_max
ie_tonic=0;
ie_pulse=0;%1e-8;
t_min=254;
t_max=263;

%% 
%%%%%%%% --------- FitzHugh-Nagumo parameters in ONCB --------- %%%%%%%%%

epsilon = 0.000002;
A = -0.01;
vz0 = -58;
vz1= -48;
vz2 = -40;%-30;%-17;%-27;%-30;
vz3 = -35;%-15;%-25;
vz4 = -61;%-63;
w0 = 0;
gamma = 30;%50;
v0 = -58.5;%-58;%-55;
exp_amp=0.3;
exp_thresh=-56;
exp_steep=0.2;
tau = 8000;%2000;%8000;%16000;%5000;
ONCB_area_factor=0.5;%0.64;%0.2;
time_end = 3000;
tau_w=2.5;
delta_tau=-.95;%-0.95;
v_tau=-55;
v_xi=-1;%-0.1;
delta_tauM=0;
v_tauM=-45;
v_xiM=-0.1;

v_shift=-1;
vz0=vz0+v_shift; vz1=vz1+v_shift;
vz2=vz2+v_shift; vz3=vz3+v_shift;
vz4=vz4+v_shift; 
v_tau=v_tau+v_shift; v0=v0+v_shift; exp_thresh=exp_thresh+v_shift;

%%
%%%%%%%% -------------------- General ---------------- %%%%%%%%%

ena = 50; %(mV)
ek = -77; %(mV)
Cm = 1e-3; %(mf/cm^2)

%%
%%%%%%%% -------------------- ONCB ------------------- %%%%%%%%%
v_init_ONCB = -60;
w_init = 0;

SA_ONCB= 439.8*1e-8; %(cm^2)
SA_ONCB= SA_ONCB*ONCB_area_factor; %(cm^2)

capacitanceONCB = Cm*SA_ONCB;


RmONCB = 12000;
epasONCB = -50;
gpasONCB=1/RmONCB;

if ONCB == 5
    gjONCB = 250*1e-12; %(S)
elseif ONCB == 7
%     gjONCB=  1250*1e-12; %(S)
    gjONCB=  750*1e-12; %(S)
elseif ONCB == 10
%     gjONCB=  1250*1e-12; %(S)
    gjONCB=  500*1e-12; %(S)
elseif ONCB == 0
    gjONCB = 0;
end

%%
%%%%%%%% -------------------- AII -------------------- %%%%%%%%%

v_init = -60;%-50;

armsection = 1; %11;

somadiam=25*1e-4;  %(cm)
somalength=25*1e-4;  %(cm)
armdiam=0.3*1e-4;  %(cm)
armlength=32*1e-4;  %(cm)
handdiam=2*1e-4;  %(cm)
handlength=2*1e-4;  %(cm)

SA_soma= 1963.5*1e-8; %(cm^2)
SA_arm=30.1593*1e-8; %(cm^2)
SA_hand=12.5664*1e-8; %(cm^2)

global_ra = 150; %(ohm-cm)

vhalfh_na = -49.5; %(mV)
vhalfm_na = -48; %(mV)
kh_na = 2; %(mV)
km_na = 5; %(mV)
mtau_na = 0.01; %(ms)
htau_na = 0.5; %(ms)

vhalfm_M = -40; %(mV) 
km_M = 4; %(mV)
mtau_M = 50;%(ms)

vhalfh_A = -40.5; %(mV)
vhalfm_A = -10; %(mV)
vhalfw_A= -45; %(mV)
kh_A = 2; %(mV)
km_A = 7; %(mV)
kw_A = 15; %(mV)
f=0.83;
mtau_A = 1; %(ms)

gbar_na_soma = 0;
gbar_na_arm = 0;
gbar_na_IS = 0.2; %(mho/cm^2) 

gbar_M_soma = 0;
gbar_M_arm = 0;
gbar_M_IS = 0.03;%0.035;%0.026;%0.03;%0.03;%0.015;%0.01;%*4/3; %(mho/cm2)
% gbar_M_IS = 0.03;


% variable g_m in IS
t_switch=1000;
t_width=1500;
gbar_M_IS_1=0.035;
gbar_M_IS_2=0.015;

figure(100)
nVAII=5;
for iVAII=1:nVAII
    VAII=-60+iVAII*5;
    VCB=-65:0.1:-15;
    gjONCB_plot=gjONCB;
    Vnull=  (tau*epsilon/capacitanceONCB) * gjONCB_plot*(VAII-VCB)...
        +A*(VCB-vz0).*(VCB-vz1).*(VCB-vz2)./(-VCB+vz3)./(VCB-vz4) - w0+exp_amp./(1+exp(-(VCB-exp_thresh)/exp_steep)); % v_ONCB
    wnull=(VCB  - v0)/gamma; %w
    subplot(1,nVAII,iVAII)
    plot(VCB,Vnull);
    hold all;
    plot(VCB,wnull);
    xlabel(' ONCB Voltage')
    ylabel(' w')
    title([' Nullclines V_{AII}=',num2str(VAII),' (mV)']);
    axis([-65 -30 -0.05 0.4])
    hold off
end
drawnow

% rd1 model has reduced A current
%gbar_A_soma = 0.004* 3/4;  %(mho/cm2) 
%gbar_A_arm = 0;
%gbar_A_IS = 0.08* 3/4;  %(mho/cm2) 

% wt has full A current
A_current_factor=1;
gbar_A_soma = 0.004*A_current_factor;  %(mho/cm2) 
gbar_A_arm = 0;
gbar_A_IS = 0.08*A_current_factor;  %(mho/cm2) 



% if numAII==1
%     Rm = 10000; %(ohm*cm^2)
%     epas = -40; %(mV)
% else
%     Rm = 40000; %(ohm*cm^2)
%     epas = -10; %(mV)
% end

%%%%%  Rd case, probably only for a network of AIIs
Rm = 40000; %(ohm*cm^2)
epas = -60;%-30; %(mV) % -60; %(mV)
gpas=1/Rm; %(S/cm^2)


gjCondR=  6000*1e-12;  %(S) : Mark's newer version (network oscillation) has 6000 pS or 4000pS, while the older version has 400 pS
gjCondIS= 1e-5*1e-12; %(S)
gjCondDC= 1e-5*1e-12; %(S)


R=zeros(armsection+1,1);
for num=1:armsection+1
    if num==1
        R(num)=(global_ra/(pi*(somadiam/2)^2))*somalength/2 + (global_ra/(pi*(armdiam/2)^2))*armlength/armsection/2;
    elseif num==armsection+1
        R(num)=(global_ra/(pi*(handdiam/2)^2))*handlength/2 + (global_ra/(pi*(armdiam/2)^2))*armlength/armsection/2;
    else
        R(num)=(global_ra/(pi*(armdiam/2)^2))*armlength/armsection;
    end
end


%%
gbar_na=zeros(numAII,armsection+2);
gbar_M=zeros(numAII,armsection+2);
gbar_A=zeros(numAII,armsection+2);
capacitance=zeros(numAII,armsection+2);
G_gap=zeros(numAII,armsection+2);
Gpas=zeros(numAII,armsection+2);

for ii=1:numAII
    
    for jj=1:armsection+2
        if jj==1
            gbar_na(ii,jj) = gbar_na_soma * SA_soma;
            gbar_M(ii,jj) = gbar_M_soma * SA_soma;
            gbar_A(ii,jj)= gbar_A_soma * SA_soma;
            capacitance(ii,jj) = Cm*SA_soma;
            G_gap(ii,jj) = gjCondR;
            Gpas(ii,jj) = gpas * SA_soma;
            
            
        elseif jj==armsection+2
            gbar_na(ii,jj) = gbar_na_IS * SA_hand;
            gbar_M(ii,jj) = gbar_M_IS * SA_hand;
            gbar_A(ii,jj) = gbar_A_IS * SA_hand;
            capacitance(ii,jj) = Cm*SA_hand;
            G_gap(ii,jj) = gjCondIS;
            Gpas(ii,jj) = gpas * SA_hand;

        else
            gbar_na(ii,jj) = gbar_na_arm * SA_arm/armsection;
            gbar_M(ii,jj) = gbar_M_arm * SA_arm/armsection;
            gbar_A(ii,jj) = gbar_A_arm * SA_arm/armsection;
            capacitance(ii,jj) = Cm*SA_arm/armsection;
            G_gap(ii,jj) = gjCondDC;
            Gpas(ii,jj) = gpas * SA_arm/armsection;
            
        end
        
    end
end


%%%%%% ------------------- Initial conditions --------------- %%%%%%%%%%
minit_na = 1/(1+exp(-(v_init-vhalfm_na)/km_na));
hinit_na = 1/(1+exp((v_init-vhalfh_na)/kh_na));
minit_M = 1/(1 + exp(-(v_init-vhalfm_M)/km_M));
minit_A = 1/(1+exp(-(v_init-vhalfm_A)/km_A));
hinit_A = f*(1/(1 + exp((v_init-vhalfh_A)/kh_A)))+(1-f);

initial = zeros(numAII, armsection+2, 9);
for ii = 1:numAII
        for jj = 1:armsection+2
            initial(ii,jj,1) = minit_na;
            initial(ii,jj,2) = hinit_na;
            initial(ii,jj,3) = minit_M;
            initial(ii,jj,4) = minit_A;
            initial(ii,jj,5) = hinit_A;
            initial(ii,jj,6) = hinit_A;
            initial(ii,jj,7) = v_init;
            initial(ii,jj,8) = v_init_ONCB;
            initial(ii,jj,9) = w_init;
            
        end
end



initial = reshape(initial,numAII*(armsection+2)*9,1);

time_start = 0;
%time_end = 800;
options = odeset('RelTol', 1e-8);
%options = odeset('RelTol', 1e-14);
[T,Y1] = ode15s(@solve_RelaxOscl_2, [time_start time_end],initial,options);
Y = zeros(numAII,armsection+2,9,length(T));
for i = 1:length(T)
    Y(1:numAII,1:armsection+2,1:9,i) = reshape(Y1(i,:),numAII,armsection+2,9);
end


plotcell = 1;
plotsection = 1; % 1 for soma, 3 for IS
v= Y(plotcell,plotsection,7,:);
v_fix=reshape(v,size(T));

plotsection = 3; % 1 for soma, 3 for IS
m_time= Y(plotcell,plotsection,3,:);
m_fix=reshape(m_time,size(T));
vs_time= Y(plotcell,plotsection,7,:);
vs_fix=reshape(vs_time,size(T));

plotsection = 1; % 1 for soma, 3 for IS
vCB= Y(plotcell,plotsection,8,:);
vCB_fix=reshape(vCB,size(T));

w= Y(plotcell,plotsection,9,:);
w_fix=reshape(w,size(T));

figure(10)
gbar_MIS_T=gbar_M_IS_1+(gbar_M_IS_2-gbar_M_IS_1)*.5*(1+tanh((T-t_switch)/t_width));
subplot(5,1,1), plot(T,v_fix), legend('V_{AII}')
%subplot(1,4,2), plot(T,vs_fix), legend('V_{AII} IS')
subplot(5,1,2), plot(T,m_fix), legend('m_{AII}')
subplot(5,1,3), plot(T,vCB_fix), legend('V_{ONCB}')
subplot(5,1,4), plot(T,w_fix), legend('w')
subplot(5,1,5), plot(T,gbar_MIS_T), legend('g_m')

figure(11)
plot(vCB_fix,w_fix)
xlabel('V_{ONCB}')
ylabel('w')
title(' ONCB Phase Plane')

end


%%%%%%%%%%%%----------------- Function called in run_RelaxOscl.m -----------------------------%%%%%%%%%%%
% Units: mF, mV,ms,S,ohm,mA, cm

function dy = solve_RelaxOscl_2(t,y1)

global epsilon A vz0 vz1 vz2 vz3 vz4 vz5 w0 gamma v0 exp_amp exp_thresh exp_steep tau_w
global delta_tau v_tau v_xi delta_tauM v_tauM v_xiM
global gjONCB ie ie_tonic ie_pulse t_min t_max
global numAII armsection R capacitance            
global gbar_na gbar_M gbar_A G_gap Gpas epas ena ek      
global vhalfh_na vhalfm_na kh_na km_na vhalfm_M km_M vhalfh_A vhalfm_A vhalfw_A kh_A km_A kw_A f
global mtau_na htau_na mtau_M mtau_A 
global tau capacitanceONCB 
global t_switch t_width gbar_M_IS_1 gbar_M_IS_2 SA_hand

y = reshape(y1,numAII,armsection+2,9);
dy = zeros(numAII,armsection+2,9);

% current injection
if (t>=t_min && t<=t_max)
    ie_t=ie_pulse;
else
    ie_t=ie_tonic;
end


gbar_MIS_t=gbar_M_IS_1+(gbar_M_IS_2-gbar_M_IS_1)*.5*(1+tanh((t-t_switch)/t_width));
gbar_M(1,armsection+2)=gbar_MIS_t * SA_hand;
%%%%%%%%%%%% ---------------------------- Differential Equations ------------------------- %%%%%%%%%%%
for i=1:numAII
    
    for j=1:armsection+2
        if (1==2)
            mtau_Mv=mtau_M*(1+delta_tauM/(1+exp(-(y(i,j,8)-v_tauM)/v_xiM)));
        else
            mtau_Mv=mtau_M;
        end
        dy(i,j,1) = (1/mtau_na)*(1/(1+exp(-(y(i,j,7)-vhalfm_na)/km_na))-y(i,j,1));  % m_na (AII)
        dy(i,j,2) = (1/htau_na)*(1/(1+exp((y(i,j,7)-vhalfh_na)/kh_na))-y(i,j,2));  %h_na (AII)
        dy(i,j,3) = (1/mtau_Mv)*(1/(1 + exp(-(y(i,j,7)-vhalfm_M)/km_M))-y(i,j,3));  %m_M (AII)
        dy(i,j,4) = (1/mtau_A)*(1/(1 + exp(-(y(i,j,7)-vhalfm_A)/km_A))-y(i,j,4));  %m_A (AII)
        dy(i,j,5) = (1/(-20/(1+exp(-(y(i,j,7)+35)/6))+25))*(f*(1/(1 + exp((y(i,j,7)-vhalfh_A)/kh_A)))+(1-f)-y(i,j,5));  %h0_A (AII)
        dy(i,j,6) = (1/min(100,((y(i,j,7)+17)^2/4+26)))*(f*(1/(1 + exp((y(i,j,7)-vhalfh_A)/kh_A)))+(1-f)-y(i,j,6));  %h1_A (AII)
        

        if j==1
            isection = 1/R(j)*(y(i,j+1,7)-y(i,j,7));
        elseif j==armsection+2
            isection = 1/R(j-1)*(y(i,j-1,7)-y(i,j,7));
        else
            isection = 1/R(j)*(y(i,j+1,7)-y(i,j,7))+1/R(j-1)*(y(i,j-1,7)-y(i,j,7));
        end
 
        if numAII == 1
            igap = 0;
        else
            if i==1
                igap = G_gap(i,j)*(y(i+1,j,7)-y(i,j,7));
            elseif i==numAII
                igap = G_gap(i,j)*(y(i-1,j,7)-y(i,j,7));
            else
                igap = G_gap(i,j)*(y(i+1,j,7)-y(i,j,7))+G_gap(i,j)*(y(i-1,j,7)-y(i,j,7));
            end
        end

        im = gbar_na(i,j)*y(i,j,1).^3.*y(i,j,2).*(y(i,j,7) - ena)+gbar_M(i,j)*y(i,j,3)*(y(i,j,7) - ek)+1/(1+exp(-(y(i,j,7)-vhalfw_A)/kw_A))*gbar_A(i,j)*y(i,j,4)*y(i,j,5)*(y(i,j,7) - ek) + (1-1/(1+exp(-(y(i,j,7)-vhalfw_A)/kw_A)))*gbar_A(i,j)*y(i,j,4)*y(i,j,6)*(y(i,j,7) - ek)+Gpas(i,j)*(y(i,j,7)-epas);
        
        if (j==1)
            i_fromONCB = gjONCB*(y(i,j,8)-y(i,j,7));
        else
            i_fromONCB=0;
        end
        
%        dy(i,j,7) = (1/capacitance(i,j))*(-im + isection + igap + i_fromONCB + ie); % V_AII
        dy(i,j,7) = (1/capacitance(i,j))*(-im + isection + igap + i_fromONCB + ie_t); % V_AII
        
        if j == 1
            coupling = 1;   %j=1 is the soma
        else
            coupling = 0;
        end
        
        if (1==1)
            dy(i,j,8) =  coupling*((1/capacitanceONCB) * gjONCB*(y(i,j,7)-y(i,j,8))...
                +(1/tau)*(1/epsilon)*(A*(y(i,j,8)-vz0).*(y(i,j,8)-vz1).*(y(i,j,8)-vz2)./(-y(i,j,8)+vz3)./(y(i,j,8)-vz4) - y(i,j,9) - w0+exp_amp/(1+exp(-(y(i,j,8)-exp_thresh)/exp_steep)))); % v_ONCB
        else
            dy(i,j,8) =  coupling*((1/capacitanceONCB) * gjONCB*(y(i,j,7)-y(i,j,8))...
                +(1/tau)*(1/epsilon)*(A*(y(i,j,8)-vz0).*(y(i,j,8)-vz1).*(y(i,j,8)-vz2).*(y(i,j,8)-vz5)./(-y(i,j,8)+vz3) - y(i,j,9) - w0+exp_amp/(1+exp(-(y(i,j,8)-exp_thresh)/exp_steep)))); % v_ONCB
        end
        if (1==1)
            dy(i,j,9) = coupling * 1/(tau*tau_w*(1+delta_tau/(1+exp(-(y(i,j,8)-v_tau)/v_xi))))*(y(i,j,8) - gamma*y(i,j,9) - v0); %w
        else
            dy(i,j,9) = coupling * (1/tau)*(y(i,j,8) - gamma*y(i,j,9) - v0); %w
        end
        
%         gjONCB
        
    end
end

dy = reshape(dy,numAII*(armsection+2)*9,1);

end

