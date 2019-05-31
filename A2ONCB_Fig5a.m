function freqVm

%%%%%%%%%%%%----------- MATLAB reconstruction of Mark's AII model on NEURON --------------%%%%%%%%%%
% Units: mF, mV,ms,S,ohm,mA, cm

clc
clear all
close all

global numcell armsection R capacitance                  % Number of cells and compartments, R and capacitance
global gbar_na gbar_M gbar_A G_gap Gpas epas ena ek      % conductances
global vhalfh_na vhalfm_na kh_na km_na                   % Na current activation/inactivation
global vhalfm_M km_M                                     % M-type K (slow) current activation
global vhalfh_A vhalfm_A vhalfw_A kh_A km_A kw_A f       % A-type K (fast) current activation/inactivation
global mtau_na htau_na mtau_M mtau_A                     % Time constants
global Iinject                                              % Injected current
global epas2

%%%%%%%%%% -------------------------- Switches -------------------------------------%%%%%%%%
IeVec = [0:-0.5:-4 -4.1:-0.05:-4.5 -5.5:-1:-6.5]*1e-9;% -4.1:-0.1:-5.5]*1e-9;
% IeVec = [-10:-0.5:-14.5 -14.6:-0.1:-16]*1e-9;
VminVec = zeros(1,length(IeVec));
VmVec = zeros(1,length(IeVec));
Vave = zeros(1,length(IeVec));

tol_steady = 0.1; %(ms)

plotswitch = 2;
lookcell = 2; % 1 for cell 1 (AII), 2 for cell 2 (CB)
numcell = 2; %Let's attach a passive cell corresponding to ON-CB very quick

plotsection =1; %1 is Soma; 13 is IS

v_init =  -45;%(mV)
v_init2 = -45;

gjCondR=  750*1e-12; %(S) %%% might want to change for CB
gjCondIS= 1e-5*1e-12; %(S)
gjCondDC= 1e-5*1e-12; %(S)

%%%%%%%%%% -------------------------- Parameters Unchanged --------------------------------%%%%%%%%
armsection = 1;%11;

somadiam=25*1e-4;  %(cm)
somalength=25*1e-4;  %(cm)
armdiam=0.3*1e-4;  %(cm)
armlength=32*1e-4;  %(cm)
handdiam=2*1e-4;  %(cm)
handlength=2*1e-4;  %(cm)

SA_soma= 1963.5*1e-8; %(cm^2)
SA_arm=30.1593*1e-8; %(cm^2)
SA_hand=12.5664*1e-8; %(cm^2)

SA_ONCB= 439.8*1e-8; %(cm^2)


global_ra = 150; %(ohm-cm)

ena = 50; %(mV)
ek = -77; %(mV)

Cm = 1e-3; %(mf/cm^2)


%%
%%% --------- Na and K currents activation/inactivation -------- %%%
vhalfh_na = -49.5; %(mV)
vhalfm_na = -48; %(mV)
kh_na = 2; %(mV)
km_na = 5; %(mV)
mtau_na = 0.01; %(ms)
htau_na = 0.5; %(ms)

vhalfm_M = -40; %(mV) % changed from -40 for test
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
gbar_na_IS = 0.2; %(mho/cm^2) changed from 0.2 to test

gbar_M_soma = 0;
gbar_M_arm = 0;
gbar_M_IS = 0.03;%*4/3;  %(mho/cm2) changed from 0.03 to test

gbar_A_soma = 0.004;%* 3/4;  %(mho/cm2) % changed from 0.004 to test
gbar_A_arm = 0;
gbar_A_IS = 0.08;%* 3/4;  %(mho/cm2) % changed from 0.08 to test


if numcell == 2
    
    gbar_na_soma2 = 0;
    gbar_na_arm2 = 0;
    gbar_na_IS2 = 0; %(mho/cm^2)
    
    gbar_M_soma2 = 0;
    gbar_M_arm2 = 0;
    gbar_M_IS2 = 0;  %(mho/cm2)
    
    gbar_A_soma2 = 0;  %(mho/cm2)
    gbar_A_arm2 = 0;
    gbar_A_IS2 = 0;  %(mho/cm2)
    
end


%%%%%%%% ----------------------------- Leak current for AII ------------------------------------ %%%%%%%
% if numcell==1
%     Rm = 10000; %(ohm*cm^2)
%     epas = -40; %(mV)
% else
%     Rm = 40000; %(ohm*cm^2)
%     epas = -10; %(mV)
% end

Rm = 40000; %(ohm*cm^2)
epas = -65; %(mV)
gpas=1/Rm; %(S/cm^2)

Rm2 = 12000;
epas2 = -35;
gpas2=1/Rm2;

%%%%%%%% --------------------------  Resistance between sections -------------------------------- %%%%%%%
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



%%%%%%%% ------------------------------ Conductances Set Up -------------------------%%%%%%%%
gbar_na=zeros(numcell,armsection+2);
gbar_M=zeros(numcell,armsection+2);
gbar_A=zeros(numcell,armsection+2);
capacitance=zeros(numcell,armsection+2);
G_gap=zeros(numcell,armsection+2);
Gpas=zeros(numcell,armsection+2);




for ii=1:numcell
    
    if ii == 2
        gbar_na_soma = gbar_na_soma2;
        gbar_na_arm = gbar_na_arm2;
        gbar_na_IS = gbar_na_IS2; %(mho/cm^2)
        
        gbar_M_soma = gbar_M_soma2;
        gbar_M_arm = gbar_M_arm2;
        gbar_M_IS = gbar_M_IS2;  %(mho/cm2)
        
        gbar_A_soma = gbar_A_soma2;  %(mho/cm2)
        gbar_A_arm = gbar_A_arm2;
        gbar_A_IS = gbar_A_IS2;  %(mho/cm2)
        
        gpas = gpas2;
        SA_soma = SA_ONCB;
    end
    
    
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
            G_gap(ii,jj) = 0;%gjCondIS;
            
            if ii == 1
                Gpas(ii,jj) = gpas * SA_hand;
            elseif ii == 2
                Gpas(ii,jj) = 0;
            end
            
            
        else
            gbar_na(ii,jj) = gbar_na_arm * SA_arm/armsection;
            gbar_M(ii,jj) = gbar_M_arm * SA_arm/armsection;
            gbar_A(ii,jj) = gbar_A_arm * SA_arm/armsection;
            capacitance(ii,jj) = Cm*SA_arm/armsection;
            G_gap(ii,jj) = 0;%gjCondDC;
            
            if ii == 1
                Gpas(ii,jj) = gpas * SA_arm/armsection;
            elseif ii == 2
                Gpas(ii,jj) = 0;
            end
            
        end
        
    end
end



%%%%%% ------------------- Initial conditions --------------- %%%%%%%%%%

%%--------------- Initial conditions for Na
minit_na = 1/(1+exp(-(v_init-vhalfm_na)/km_na));
hinit_na = 1/(1+exp((v_init-vhalfh_na)/kh_na));

%%--------------- Initial conditions for M-type K
minit_M = 1/(1 + exp(-(v_init-vhalfm_M)/km_M));

%%--------------- Initial conditions for A-type K
minit_A = 1/(1+exp(-(v_init-vhalfm_A)/km_A));
hinit_A = f*(1/(1 + exp((v_init-vhalfh_A)/kh_A)))+(1-f);

%%%%%% ------------------- Initial conditions for cell 2 --------------- %%%%%%%%%%
if numcell == 2
    minit_na2 = 1/(1+exp(-(v_init2-vhalfm_na)/km_na));
    hinit_na2 = 1/(1+exp((v_init2-vhalfh_na)/kh_na));
    
    minit_M2 = 1/(1 + exp(-(v_init2-vhalfm_M)/km_M));
    
    minit_A2 = 1/(1+exp(-(v_init2-vhalfm_A)/km_A));
    hinit_A2 = f*(1/(1 + exp((v_init2-vhalfh_A)/kh_A)))+(1-f);
end

%%%%%% ------------------- Initial conditions Set Up  ------------------ %%%%%%%%%%
initial = zeros(numcell, armsection+2, 7);
for ii = 1:numcell
    if ii==1
        for jj = 1:armsection+2
            initial(ii,jj,1) = minit_na;
            initial(ii,jj,2) = hinit_na;
            initial(ii,jj,3) = minit_M;
            initial(ii,jj,4) = minit_A;
            initial(ii,jj,5) = hinit_A;
            initial(ii,jj,6) = hinit_A;
            initial(ii,jj,7) = v_init;
        end
    elseif ii==2
        for jj = 1:armsection+2
            initial(ii,jj,1) = minit_na2;
            initial(ii,jj,2) = hinit_na2;
            initial(ii,jj,3) = minit_M2;
            initial(ii,jj,4) = minit_A2;
            initial(ii,jj,5) = hinit_A2;
            initial(ii,jj,6) = hinit_A2;
            initial(ii,jj,7) = v_init2;
        end
    end
end
initial = reshape(initial,numcell*(armsection+2)*7,1);


vmaxmat = zeros(1,length(IeVec));%[];
vminmat = zeros(1,length(IeVec));%[];
frequencymat = zeros(1,length(IeVec));%[];



for ii = 1:length(IeVec)
    
    Iinject = IeVec(ii);
    
    time_start = 0;
    time_end = 1500;
    cuttime = 500;
        
        %%%%%%%%%%% ------------------- Solving Diff Eqs ------------------------ %%%%%%%%%%
        
        options = odeset('RelTol', 1e-8);
        [T,Y1] = ode15s(@solve_FreqVm, [time_start time_end], initial,options);
        Y = zeros(numcell,armsection+2,7,length(T));
        for i = 1:length(T)
            Y(1:numcell,1:armsection+2,1:7,i) = reshape(Y1(i,:),numcell, armsection+2,7);
        end
        initial_new = reshape(Y(:,:,:,end),numcell*(armsection+2)*7,1);
        initial = initial_new;
        
        plotcell =1;
        v= Y(plotcell,plotsection,7,:);
        v_fix=reshape(v,size(T));
        
        VminVec(ii) = min(v_fix);
        Vave(ii) = (min(v_fix)+max(v_fix))/2;

        
        %%%%% -------------- If there are more than one cell -------------- %%%%%%%
        if numcell == 2
            plotcell2 = 2;
            v2= Y(plotcell2,plotsection,7,:);
            v2_fix=reshape(v2,size(T));
        end
        
        %%% ----------------- Cutting of the transient ------------------%%%
        savenow = [];
        for i = 1:length(T)
            if T(i) >= cuttime
                savenow = [savenow i];
            else
                savenow = savenow;
            end
            start = min(savenow);
            
        end
        
        if lookcell == 1
            v = v_fix(start:end);
        elseif lookcell == 2
            v = v2_fix(start:end);
        end
        
        t = T(start:end);
        
        
        %%% ----------------- Finding local max/minima ----------------------%%%
        if lookcell == 1
            [vmax,imax,vmin,imin] = extrema(v_fix(start:end));
        elseif lookcell == 2
            [vmax,imax,vmin,imin] = extrema(v2_fix(start:end));
        end
        
        
        imax = sort(imax);
        imin = sort(imin);
        
        
        
        %% No spiking nor bursting
        if length(vmax)<=2 || length(vmax)<=2
            
                tmax = T(imax);
                Vmax = vmax;
                tmin = T(imin);
                Vmin = vmin;
                if length(vmax)<1 || length(vmax)<1
                    freq = 0; %(Hz)
                    t_cross = [];
                else
                    freq = 0; %(Hz)
                    t_cross = [];
                end

        else
            
            %% Quadratic interpolation
            
            tmax = zeros(1,length(imax));
            Vmax = zeros(1,length(imax));
            iMax = zeros(1,length(imax));
            
            for index = 1:length(imax)
                iMax = imax(index);
                x = [t(iMax-1) t(iMax) t(iMax+1)];
                y = [v(iMax-1) v(iMax) v(iMax+1)];
                
                [p,S,mu] = polyfit(x,y,2);
                A_new = p(1); B_new = p(2); C_new = p(3); mu1 = mu(1); mu2 = mu(2);
                c=A_new*(mu1)^2/(mu2)^2-(B_new*(mu1/mu2))+C_new;
                b = (B_new*mu2 - 2*A_new*mu1)/(mu2)^2;
                a = A_new/(mu2)^2;
                
                
                tmax(index) = -b/(2*a);
                Vmax(index) = c-b^2/(4*a);
            end
            
            tmin = zeros(1,length(imin));
            Vmin = zeros(1,length(imin));
            iMin = zeros(1,length(imin));
            
            for index = 1:length(imin)
                iMin = imin(index);
                x = [t(iMin-1) t(iMin) t(iMin+1)];
                y = [v(iMin-1) v(iMin) v(iMin+1)];
                
                [p,S,mu] = polyfit(x,y,2);
                A_new = p(1); B_new = p(2); C_new = p(3); mu1 = mu(1); mu2 = mu(2);
                c=A_new*(mu1)^2/(mu2)^2-(B_new*(mu1/mu2))+C_new;
                b = (B_new*mu2 - 2*A_new*mu1)/(mu2)^2;
                a = A_new/(mu2)^2;
                
                tmin(index) = -b/(2*a);
                Vmin(index) = c-b^2/(4*a);
            end
            
            %%
            crossline = (max(Vmax)+min(Vmin))/2;
            count1 = 1;
            count2 = 1;
            Vmin_re = [];
            while count1 <= length(tmin)
                if Vmin(count1)<crossline
                    Vmin_re(count2) = Vmin(count1);
                    count2=count2+1;
                end
                count1=count1+1;
            end
            
            count1 = 1;
            count2 = 1;
            Vmax_re = [];
            while count1 <= length(tmax)
                if Vmax(count1)>crossline
                    Vmax_re(count2) = Vmax(count1);
                    count2=count2+1;
                end
                count1=count1+1;
            end
            
            Vmax = Vmax_re;
            Vmin = Vmin_re;
            
            %%
            t_cross = [];
            v_cross = (max(vmax)+min(vmin))/2;
            for iter = 1:length(v)-1
                if (v(iter)-v_cross)<=0 && (v(iter+1)-v_cross)>=0
                    t_cross = [t_cross ((v_cross-v(iter))*(t(iter+1)-t(iter))/(v(iter+1)-v(iter)))+t(iter)];
                else
                    t_cross = t_cross;
                end
            end
            
            V_cross = ones(1,length(t_cross))*v_cross;
            
            if length(t_cross)<=2
                
                    
                    if length(t_cross)<=1
                        freq = 0; %(Hz)
                    else
                        period = (t_cross(2)-t_cross(1))*1e-3;
                        freq=1/period;
                    end
     
                
            else
                

                dtcross = [];
                for step = 1:length(t_cross)-1
                    dtcross_temp = t_cross(step+1)-t_cross(step);
                    if dtcross_temp>(1/30)*1000
                        dtcross = [dtcross dtcross_temp];
                    else
                        dtcross = dtcross;
                    end
                end
                
                vari = zeros(1,length(dtcross)-1);
                for step = 1:length(dtcross)-1
                    vari(step) = abs(dtcross(step+1)-dtcross(step));
                    if vari(step) < tol_steady
                        period = dtcross(step+1)*1e-3; %(s)
                        freq = 1/period; %(Hz)
                        break
                    end
                    
                end
                
                
                
                %% Increase the simulataion time if :
                % 1) it does not go to a steady state
                % 2) it does not have more than one period in it
                % Need a limit
                
                if min(vari)>=tol_steady
                    
                        period = dtcross(end)*1e-3; %(s)
                        freq = 1/period; %(Hz)
         
                end
                
            end
            
        end

    
        
    if length(Vmin)<1
        VmVec(ii) = VminVec(ii);
    else
        VmVec(ii) = min(Vmin);
    end
    
    if freq == 0
     VminVec(ii) = v_fix(end);
    end
    
    frequencymat(ii) = freq;
    
    %%%%%%%%%%------------------------------------  Decide on Plots -------------------------------------%%%%%%%%%%
    if plotswitch == 1
        
        figure(21)
        clf
        %             hold on
        plot(T,v_fix,'r')%,T,v2_fix,'b')
%         hold on,plot(tmax,Vmax,'.')
        if length(t_cross)>=1
            hold on,plot(t_cross,V_cross,'g.')
        end
        xlabel('time (ms)'), ylabel('voltage (mV)')
        
      
        pause
        
    end
    
    Vave(ii)
    VmVec(ii)
    VminVec(ii)
    freq
    Iinject
end


if plotswitch == 2
    
    figure()
    plot(VminVec,frequencymat,'.')
    xlabel('V_{min} (mV)'), ylabel('Frequency (Hz)')
    title('Frequency')

%     plot(VmVec,frequencymat,'.')
%     xlabel('V_{m}'), ylabel('Frequency (Hz)')

    figure()
    plot(Vave,frequencymat,'.')
    xlabel('V_{Ave} (mV)'), ylabel('Frequency (Hz)')
    
    figure()
    plot(IeVec.*1e9,frequencymat,'.')
    xlabel('I_{inject} (pA)'), ylabel('Frequency (Hz)')
%     figure()
%     ord = sort(VminVec,'descend');
%     plot(ord,frequencymat,'.')


    
    
    %     hgsave('Vm_func.fig');
    
end

end

%%%%%%%%%%%%----------- MATLAB reconstruction of Mark's AII model on NEURON --------------%%%%%%%%%%
% Units: mF, mV,ms,S,ohm,mA, cm


function dy = solve_FreqVm(t,y1)

global numcell armsection R capacitance                  % Number of cells and compartments, R and capacitance
global gbar_na gbar_M gbar_A G_gap Gpas epas ena ek      % conductances

global vhalfh_na vhalfm_na kh_na km_na                   % Na current activation/inactivation
global vhalfm_M km_M                                     % M-type K (slow) current activation
global vhalfh_A vhalfm_A vhalfw_A kh_A km_A kw_A f       % A-type K (fast) current activation/inactivation
global mtau_na htau_na mtau_M mtau_A                     % Time constants
global Iinject                                            % Injected current
global epas2

y = reshape(y1,numcell,armsection+2,7);
dy = zeros(numcell,armsection+2,7);


%%%%%%%%%%%% ---------------------------- Differential Equations ------------------------- %%%%%%%%%%%
for i=1:numcell

    for j=1:armsection+2
        
%         if t>=cuttime && 
        if i ==1 && j == 1 % i == 17
            ie =  Iinject;
        else
            ie = 0;
        end
        
        
        dy(i,j,1) = (1/mtau_na)*(1/(1+exp(-(y(i,j,7)-vhalfm_na)/km_na))-y(i,j,1));  % m_na
        dy(i,j,2) = (1/htau_na)*(1/(1+exp((y(i,j,7)-vhalfh_na)/kh_na))-y(i,j,2));  %h_na
        dy(i,j,3) = (1/mtau_M)*(1/(1 + exp(-(y(i,j,7)-vhalfm_M)/km_M))-y(i,j,3));  %m_M
        dy(i,j,4) = (1/mtau_A)*(1/(1 + exp(-(y(i,j,7)-vhalfm_A)/km_A))-y(i,j,4));  %m_A
        dy(i,j,5) = (1/(-20/(1+exp(-(y(i,j,7)+35)/6))+25))*(f*(1/(1 + exp((y(i,j,7)-vhalfh_A)/kh_A)))+(1-f)-y(i,j,5));  %h0_A
        dy(i,j,6) = (1/min(100,((y(i,j,7)+17)^2/4+26)))*(f*(1/(1 + exp((y(i,j,7)-vhalfh_A)/kh_A)))+(1-f)-y(i,j,6));  %h1_A
        
        if j==1
          
            isection = 1/R(j)*(y(i,j+1,7)-y(i,j,7));
            
            
        elseif j==armsection+2
            
            
            isection = 1/R(j-1)*(y(i,j-1,7)-y(i,j,7));
            
        else
                    
            
            isection = 1/R(j)*(y(i,j+1,7)-y(i,j,7))+1/R(j-1)*(y(i,j-1,7)-y(i,j,7));
            
        end
        
  
        if numcell == 1
            igap = 0;
        else
            if i==1
                igap = G_gap(i,j)*(y(i+1,j,7)-y(i,j,7));
            elseif i==numcell
                igap = G_gap(i,j)*(y(i-1,j,7)-y(i,j,7));
            else
                igap = G_gap(i,j)*(y(i+1,j,7)-y(i,j,7))+G_gap(i,j)*(y(i-1,j,7)-y(i,j,7));
            end
        end
        
        if i==2
            im = gbar_na(i,j)*y(i,j,1).^3.*y(i,j,2).*(y(i,j,7) - ena)+gbar_M(i,j)*y(i,j,3)*(y(i,j,7) - ek)+1/(1+exp(-(y(i,j,7)-vhalfw_A)/kw_A))*gbar_A(i,j)*y(i,j,4)*y(i,j,5)*(y(i,j,7) - ek) + (1-1/(1+exp(-(y(i,j,7)-vhalfw_A)/kw_A)))*gbar_A(i,j)*y(i,j,4)*y(i,j,6)*(y(i,j,7) - ek)+Gpas(i,j)*(y(i,j,7)-epas2);
            
        elseif i == 1
            
            im = gbar_na(i,j)*y(i,j,1).^3.*y(i,j,2).*(y(i,j,7) - ena)+gbar_M(i,j)*y(i,j,3)*(y(i,j,7) - ek)+1/(1+exp(-(y(i,j,7)-vhalfw_A)/kw_A))*gbar_A(i,j)*y(i,j,4)*y(i,j,5)*(y(i,j,7) - ek) + (1-1/(1+exp(-(y(i,j,7)-vhalfw_A)/kw_A)))*gbar_A(i,j)*y(i,j,4)*y(i,j,6)*(y(i,j,7) - ek)+Gpas(i,j)*(y(i,j,7)-epas);
            
        end
        
        dy(i,j,7) = (1/capacitance(i,j))*(-im + isection + igap +ie);

        
    end
end

dy = reshape(dy,numcell*(armsection+2)*7,1);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xmax,imax,xmin,imin] = extrema(x)
%EXTREMA   Gets the global extrema points from a time series.
%   [XMAX,IMAX,XMIN,IMIN] = EXTREMA(X) returns the global minima and maxima 
%   points of the vector X ignoring NaN's, where
%    XMAX - maxima points in descending order
%    IMAX - indexes of the XMAX
%    XMIN - minima points in descending order
%    IMIN - indexes of the XMIN
%
%   DEFINITION (from http://en.wikipedia.org/wiki/Maxima_and_minima):
%   In mathematics, maxima and minima, also known as extrema, are points in
%   the domain of a function at which the function takes a largest value
%   (maximum) or smallest value (minimum), either within a given
%   neighbourhood (local extrema) or on the function domain in its entirety
%   (global extrema).
%
%   Example:
%      x = 2*pi*linspace(-1,1);
%      y = cos(x) - 0.5 + 0.5*rand(size(x)); y(40:45) = 1.85; y(50:53)=NaN;
%      [ymax,imax,ymin,imin] = extrema(y);
%      plot(x,y,x(imax),ymax,'g.',x(imin),ymin,'r.')
%
%   See also EXTREMA2, MAX, MIN

%   Written by
%   Lic. on Physics Carlos Adrián Vargas Aguilera
%   Physical Oceanography MS candidate
%   UNIVERSIDAD DE GUADALAJARA 
%   Mexico, 2004
%
%   nubeobscura@hotmail.com

% From       : http://www.mathworks.com/matlabcentral/fileexchange
% File ID    : 12275
% Submited at: 2006-09-14
% 2006-11-11 : English translation from spanish. 
% 2006-11-17 : Accept NaN's.
% 2007-04-09 : Change name to MAXIMA, and definition added.


xmax = [];
imax = [];
xmin = [];
imin = [];

% Vector input?
Nt = numel(x);
if Nt ~= length(x)
 error('Entry must be a vector.')
end

% NaN's:
inan = find(isnan(x));
indx = 1:Nt;
if ~isempty(inan)
 indx(inan) = [];
 x(inan) = [];
 Nt = length(x);
end

% Difference between subsequent elements:
dx = diff(x);

% Is an horizontal line?
if ~any(dx)
 return
end

% Flat peaks? Put the middle element:
a = find(dx~=0);              % Indexes where x changes
lm = find(diff(a)~=1) + 1;    % Indexes where a do not changes
d = a(lm) - a(lm-1);          % Number of elements in the flat peak
a(lm) = a(lm) - floor(d/2);   % Save middle elements
a(end+1) = Nt;

% Peaks?
xa  = x(a);             % Serie without flat peaks
b = (diff(xa) > 0);     % 1  =>  positive slopes (minima begin)  
                        % 0  =>  negative slopes (maxima begin)
xb  = diff(b);          % -1 =>  maxima indexes (but one) 
                        % +1 =>  minima indexes (but one)
imax = find(xb == -1) + 1; % maxima indexes
imin = find(xb == +1) + 1; % minima indexes
imax = a(imax);
imin = a(imin);

nmaxi = length(imax);
nmini = length(imin);                

%%%% ----- Commented Out on 052012 by Hannah ------------------%%%%%
% Maximum or minumim on a flat peak at the ends?
if (nmaxi==0) && (nmini==0)
 if x(1) > x(Nt)
  xmax = x(1);
  imax = indx(1);
  xmin = x(Nt);
  imin = indx(Nt);
 elseif x(1) < x(Nt)
  xmax = x(Nt);
  imax = indx(Nt);
  xmin = x(1);
  imin = indx(1);
 end
 return
end

%% Maximum or minumim at the ends?
% if (nmaxi==0) 
%  imax(1:2) = [1 Nt];
% elseif (nmini==0)
%  imin(1:2) = [1 Nt];
% else
%  if imax(1) < imin(1)
%   imin(2:nmini+1) = imin;
%   imin(1) = 1;
%  else
%   imax(2:nmaxi+1) = imax;
%   imax(1) = 1;
%  end
%  if imax(end) > imin(end)
%   imin(end+1) = Nt;
%  else
%   imax(end+1) = Nt;
%  end
% end
%%%%% -----------------------------------------------------------%%%%%
xmax = x(imax);
xmin = x(imin);

% NaN's:
if ~isempty(inan)
 imax = indx(imax);
 imin = indx(imin);
end

% Same size as x:
imax = reshape(imax,size(xmax));
imin = reshape(imin,size(xmin));

% Descending order:
[temp,inmax] = sort(-xmax); clear temp
xmax = xmax(inmax);
imax = imax(inmax);
[xmin,inmin] = sort(xmin);
imin = imin(inmin);


% Carlos Adrián Vargas Aguilera. nubeobscura@hotmail.com

end

