function [dY, dati] = AVNMacromodel_new(time, Y,var)
%AVN model
% time (millisecond)
% current in pA/pF
% concentration in mM

%% Constants

F  = 96485000000000000; % femtocoulomb_per_millimole
F1 = F*(10^-12);
R  = 8314400000000000; % attojoule_per_millimole_kelvin
R2 = R*(10^-12);
T = 310; % K

V  = Y(1); %mV
%vffrt = V*F1*F1/(R2*T);
%vfrt  = V*F1/(R2*T);
%% Cell dimension

Cm      =  22;     %pF
Rcell   = 0.00005; %dm
Cm_surf = 1;       %uF/cm^2
Lcell = (1/(2*pi*Rcell))*(Cm/(10^8*Cm_surf)-(2*pi*(Rcell)^2)); %dm % 65 um

ds = 0.01; %dimensionless

%litre
V_Cell = 3.1416*Rcell^2*Lcell;
V_sub  = ds*V_Cell;
V_i    = 0.46*V_Cell-V_sub;
V_JSR  = 0.0012*V_Cell;
V_NSR  = 0.0116*V_Cell;

%% Ionic concentrations (mM)
% Extracellular concentrations
Nao = 140;
Ko  = 5.4;
Cao = 2;
% Intracellular concentrations
Nai = Y(28);
Ki  = Y(29);
Cai = Y(30);
Ca_sub = Y(31);
CaJSR  = Y(33);
CaNSR  = Y(32);

%% Nernst potential (mV)
E_Na  = (R*T)/F*log(Nao/Nai);
E_Ca  = 0.5*(R*T)/F*log(Cao/Ca_sub);
% PhiCaL=4.0.*vffrt*(1*Ca_sub.*exp(2.0*vfrt)-0.341*Cao)./(exp(2.0*vfrt)-1.0);
%E_Ca_exp = 42;
E_K   = (R*T)/F*log(Ko/Ki);
PkNa  = 0.1; 
E_K_s = (R*T)/F*log((Ko+PkNa*Nao)/(Ki+PkNa*Nai));

% block factor
inas_block = 1;
ina_block = 1;
icaD_block = 1;
ical_block = 1;
icat_block = 1;
if_block = 1; 
ikr_block = 1;
ibCa_block = 1; 
ibNa_block = 1;
ibK_block = 1;
jrel_block = 1;
pup_block = 1;

% scaling factor
scaling_INa = 2.2; 
scaling_INas = 1.5*0.17; 
scaling_IKr = 0.3; 
scaling_ICaL = 1; 
scaling_IKs = 1*80;
scaling_ik1 = 1;
scaling_ist = .2;
block_fourAP = 1;
scaling_Ito = 1.8*block_fourAP; 
scaling_Isus = 5*block_fourAP;  


%**********************************************************************************************
%  Macrophage
%*********************************************************************************************
Choix_Type_Phi=var(2);
Cm_Phi=var(14); %pF 27.9
G_gap=var(1); %0,1,2,8,nS
K_e_Phi=var(3);
%ik1_EK=-79.362; % mV; 

  %gb_Phi=0.16/Cm_Phi;
  EK_Phi=-79.362; %mV;
  F_Phi=96487;
  R_Phi=8314;
  T_Phi=297; 

  %gb_Phi=1.6e-7/Cm_Phi;

%% temperature factor
q10 = 0.2378;

%% Stimulation
% 
% i_stim_Amplitude 		= 0.3*1e3;   %nanoampere (in stim_mode)
% i_stim_End 				= 1000000.0;   % msecond (in stim_mode)
% i_stim_PulseDuration	= 3;   % msecond (in stim_mode)
% i_stim_Start 			= 0.0;   % msecond (in stim_mode)
% i_stim_frequency 		= bpm;   % beat per_minut (in stim_mode)
% stim_flag 				= 0;   % dimensionless (in stim_mode)
% i_stim_Period 			= 60.0*1e3/i_stim_frequency; %ms

i_stim_Amplitude 		=var(4);   %nanoampere (in stim_mode)
 i_stim_Amplitude_Macro=var(13);
i_stim_End 				= var(5);   % msecond (in stim_mode)
i_stim_PulseDuration	= var(6);   % msecond (in stim_mode)
i_stim_Start 			= var(7);   % msecond (in stim_mode)
%i_stim_frequency 		= var(8);   % beat per_minut (in stim_mode)
stim_flag 				=var(9);   % dimensionless (in stim_mode)
stim_flagM=var(10);
i_stim_Period 			= var(11); %ms
number_of_macrophages=var(12);


if ((time >= i_stim_Start) && (time <= i_stim_End) && (time-i_stim_Start-floor((time-i_stim_Start)/i_stim_Period)*i_stim_Period <= i_stim_PulseDuration))
      i_stim =stim_flag*i_stim_Amplitude/Cm ;
       i_diff_Phi=  stim_flagM* i_stim_Amplitude_Macro/Cm_Phi;
else
    i_stim = 0.0;
      i_diff_Phi=0.0;
end


%****************************************************************************************
  % Macrophage
  %***************************************************************************************
    

    I_Gap_Phi=(G_gap/Cm_Phi)*(Y(51)-V);
    I_Gap_AV=(G_gap/Cm)*(V-Y(51));
if  Choix_Type_Phi==0
    %**********************************************************************************************
    % Type 0
    %**********************************************************************************************
 gb_Phi=(1.6e-7/Cm_Phi)*1000000;
 GKir_Phi=(4.2e-7/Cm_Phi)*1000000;
    a_Kir=0.91;
    b_Kir=1.09;
  Eb_Phi=-10; 
    % <component name="ib">
  
    Ib_Phi=gb_Phi*(Y(51)-Eb_Phi);
     frt=F_Phi/(R_Phi*T_Phi);
    %<component name="ik1">
    IK1_Phi=GKir_Phi*sqrt(K_e_Phi)*(Y(51)-EK_Phi)/(a_Kir+exp(b_Kir*((Y(51)-EK_Phi)*frt)));

    %T=(275+(22 kelvin));
    dY(42,1)=0;
    dY(43,1)=0;
    dY(44,1)=0;
    dY(45,1)=0;
    dY(46,1)=0;
    % 
    dY(48,1)=0;
    dY(47,1)=0;
    % 
    dY(49,1)=0;
    dY(50,1)=0;

    % <component name="cell">
   dY(51,1)=-(Ib_Phi+IK1_Phi+I_Gap_Phi)+i_diff_Phi;
    Ishk_Phi=0* Y(51);
     IKur_Phi=0*Y(51);
elseif Choix_Type_Phi==1

    %******************************************************************************************************
    %Type 1
    %******************************************************************************************************
     CO=Y(42); % 
     C1=Y(43); % 
     C2=Y(44); %
     C3=Y(45); % 
     C4=Y(46); % 
     Oshaker=Y(47); %
      I=Y(48); % 
    %Y(49) % iKur.R
    %Y(50) % S
    gb_Phi=(2.2e-7/Cm_Phi)*1000000;

   GKir_Phi=(7.3e-7/Cm_Phi)*1000000;
   
   a_Kir=0.8;
   b_Kir=0.7;
   A_Phi=127.7;
   B_Phi=16.5;
   a_Phi=3.96;
   b_Phi=0.1;
   c_Phi=0.0029;
   d_Phi=0.0008;
   m_Phi=11.1;
   n_Phi=26.0;
   p_Phi=1205.3;
   q_Phi=3377.2;

    %      a_Kir=0.8*1000;
    % b_Kir=0.7*1000;
    % A_Phi=127.7*1000;
    % B_Phi=16.5*1000;
    % a_Phi=3.96*1000;
    % b_Phi=0.1*1000;
    % c_Phi=0.0029*1000;
    % d_Phi=0.0008*1000;
    % m_Phi=11.1;
    % n_Phi=26.0;
    % p_Phi=1205.3;
    % q_Phi=3377.2;
  gshaker_Phi=(0.015/Cm_Phi);
    nchannels=74;
     Eb_Phi=0;
    
    % <component name="ib">
    Ib_Phi=(gb_Phi*( Y(51)-Eb_Phi));
    frt=F_Phi/(R_Phi*T_Phi);
    %<component name="ik1">
    IK1_Phi=GKir_Phi*sqrt(K_e_Phi)*(Y(51)-EK_Phi)/(a_Kir+exp(b_Kir*((Y(51)-EK_Phi)*frt)));
    %T=(275+(22 kelvin));

    %<component name="ishk">
    alpha_Phi=(a_Phi*exp(Y(51)/m_Phi));
    beta_Phi=(b_Phi*exp((-1)*Y(51)/n_Phi));
    mu_Phi=(c_Phi*exp(Y(51)/p_Phi));
    phi_Phi=(d_Phi*exp((-1)*Y(51)/q_Phi));



    dY(42,1)=((-1)*4*alpha_Phi*CO+beta_Phi*C1);
    dY(43,1)=(4*alpha_Phi*CO-(3*alpha_Phi+beta_Phi)*C1+2*beta_Phi*C2);
    dY(44,1)=(3*alpha_Phi*C1-(2*alpha_Phi+2*beta_Phi)*C2+3*beta_Phi*C3);
    dY(45,1)=2*alpha_Phi*C2-(alpha_Phi+3*beta_Phi)*C3+4*beta_Phi*C4;
    dY(46,1)=(alpha_Phi*C3-(A_Phi+4*beta_Phi)*C4+B_Phi*Oshaker);
    ishk_EK=-79.362;% mV
    dY(48,1)=(mu_Phi*Oshaker-phi_Phi*I);
    dY(47,1)=(A_Phi*C4-(mu_Phi+B_Phi)*Oshaker+phi_Phi*I);
    Ishk_Phi=(gshaker_Phi*nchannels*Oshaker*( Y(51)-ishk_EK));
       dY(49,1)=0;
    dY(50,1)=0;
 
    dY(51,1)=-(Ishk_Phi+Ib_Phi+IK1_Phi+I_Gap_Phi)+i_diff_Phi;
      IKur_Phi=0*Y(51);
   elseif Choix_Type_Phi==2
    %*******************************************************************************************************
    %Type2
    %*********************************************************************************************************
     GKur_Phi=(3e-7/Cm_Phi)*1000000;
   gb_Phi=(3.5e-7/Cm_Phi)*1000000;

 GKir_Phi=(1.3e-6/Cm_Phi)*1000000;

    a_Kir=0.72;
    b_Kir=0.96;
    aKur_Phi=86;
    bKur_Phi=18.8;
    cKur_Phi=1.2;
    dKur_Phi=1;
   Eb_Phi=8;
    %<component name="ib">
    Ib_Phi=(gb_Phi*( Y(51)-Eb_Phi));

    % <component name="iKur">
    iKur_EK=-79.362; % mV
    IKrinf=(1/(1+exp(( Y(51)+6)/aKur_Phi)));
    IKrtau=(.0066/(1+exp(( Y(51)+5)/cKur_Phi))+3.6e-4);
    IKsinf=(1/(1+exp(( Y(51)+7.5)/bKur_Phi)));
    IKstau=(.43/(1+exp(( Y(51)+60)/dKur_Phi))+2.2);
    IKur_Phi=(GKur_Phi*Y(49)*Y(50)*( Y(51)-iKur_EK));
    dY(42,1)=0;
    dY(43,1)=0;
    dY(44,1)=0;
    dY(45,1)=0;
    dY(46,1)=0;
 
    dY(48,1)=0;
    dY(47,1)=0;

    dY(49,1)=((IKrinf-Y(49))/IKrtau);
    dY(50,1)=((IKsinf-Y(50))/IKstau);
    



     frt=F_Phi/(R_Phi*T_Phi);
    %<component name="ik1">
    IK1_Phi=GKir_Phi*sqrt(K_e_Phi)*(Y(51)-EK_Phi)/(a_Kir+exp(b_Kir*((Y(51)-EK_Phi)*frt)));

    %<component name="cell">
    dY(51,1)=-(IKur_Phi+Ib_Phi+IK1_Phi+I_Gap_Phi)+i_diff_Phi;
   Ishk_Phi= 0*Y(51);

elseif Choix_Type_Phi==12
 %*******************************************************************************************************
    %Type 12
    %*********************************************************************************************************
 GKur_Phi=(3e-7/Cm_Phi)*1000000;
   gb_Phi=(3.5e-7/Cm_Phi)*1000000;

 GKir_Phi=(1.3e-6/Cm_Phi)*1000000;
    CO=Y(42); % 
     C1=Y(43); % 
     C2=Y(44); %
     C3=Y(45); % 
     C4=Y(46); % 
     Oshaker=Y(47); %
      I=Y(48); % 

      
    a_Kir=0.72;
    b_Kir=0.96;
    aKur_Phi=86;
    bKur_Phi=18.8;
    cKur_Phi=1.2;
    dKur_Phi=1;

     A_Phi=127.7;
    B_Phi=16.5;
    a_Phi=3.96;
    b_Phi=0.1;
    c_Phi=0.0029;
    d_Phi=0.0008;
    m_Phi=11.1;
    n_Phi=26.0;
    p_Phi=1205.3;
    q_Phi=3377.2;
    nchannels=74;
     gshaker_Phi=0.015/Cm_Phi;

   Eb_Phi=8;
    %<component name="ib">
    Ib_Phi=(gb_Phi*( Y(51)-Eb_Phi));

    % <component name="iKur">
    iKur_EK=-79.362; % mV
    IKrinf=(1/(1+exp(( Y(51)+6)/aKur_Phi)));
    IKrtau=(.0066/(1+exp(( Y(51)+5)/cKur_Phi))+3.6e-4);
    IKsinf=(1/(1+exp(( Y(51)+7.5)/bKur_Phi)));
    IKstau=(.43/(1+exp(( Y(51)+60)/dKur_Phi))+2.2);
    IKur_Phi=(GKur_Phi*Y(49)*Y(50)*( Y(51)-iKur_EK));
  

     dY(49,1)=((IKrinf-Y(49))/IKrtau);
    dY(50,1)=((IKsinf-Y(50))/IKstau);
    



     frt=F_Phi/(R_Phi*T_Phi);
    %<component name="ik1">
    IK1_Phi=GKir_Phi*sqrt(K_e_Phi)*(Y(51)-EK_Phi)/(a_Kir+exp(b_Kir*((Y(51)-EK_Phi)*frt)));

   %<component name="ishk">
    alpha_Phi=(a_Phi*exp(Y(51)/m_Phi));
    beta_Phi=(b_Phi*exp((-1)*Y(51)/n_Phi));
    mu_Phi=(c_Phi*exp(Y(51)/p_Phi));
    phi_Phi=(d_Phi*exp((-1)*Y(51)/q_Phi));



     dY(42,1)=((-1)*4*alpha_Phi*CO+beta_Phi*C1);
    dY(43,1)=(4*alpha_Phi*CO-(3*alpha_Phi+beta_Phi)*C1+2*beta_Phi*C2);
    dY(44,1)=(3*alpha_Phi*C1-(2*alpha_Phi+2*beta_Phi)*C2+3*beta_Phi*C3);
    dY(45,1)=2*alpha_Phi*C2-(alpha_Phi+3*beta_Phi)*C3+4*beta_Phi*C4;
    dY(46,1)=(alpha_Phi*C3-(A_Phi+4*beta_Phi)*C4+B_Phi*Oshaker);
    ishk_EK=-79.362;% mV
    dY(48,1)=(mu_Phi*Oshaker-phi_Phi*I);
    dY(47,1)=(A_Phi*C4-(mu_Phi+B_Phi)*Oshaker+phi_Phi*I);
    Ishk_Phi=(gshaker_Phi*nchannels*Oshaker*( Y(51)-ishk_EK));
  




    %<component name="cell">
    dY(51,1)=-(IKur_Phi+Ib_Phi+IK1_Phi+ Ishk_Phi+I_Gap_Phi)+i_diff_Phi;
  

 end

%
%
%**********************************************************************************************


%% INas

ms  = Y(5);
hs1 = Y(6);
hs2 = Y(7);
jhs = Y(41);

g_Nas   = (2.5e-06/Cm*1000)*scaling_INas; %nS/pF

INas_ac = 38; %mV reform
sNas    = 6; %reform
INas_In = -56; %mV

ms_inf  = (1.0/(1+exp((V+INas_ac)/(-sNas))));
tau_ms  = q10*((0.6247/(0.832*exp(-0.335*(V+56.7))+0.627*exp(0.082*(V+65.01))))+0.04); %ms
dY(5,1) = (ms_inf-ms)/tau_ms;

hs1_inf = 1/(1+exp((V-INas_In)/3));
tau_hs1=q10*0.113*(1/(13475.066*exp((V-59.398)/15.645))+1/(1.113+0.044*exp(-(V-86.768)/8.059)))^(-1); % KHARCHE
dY(6,1) = (hs1_inf-hs1)/tau_hs1;

hs2_inf = hs1_inf;
par1 = 0.125;
par2 = 140557.232; par3 = 59.455; par4 = 12; 
par5 = 2.471; par6 = 0.767; par7 = 68.931; par8 = 18.237;
tau_hs2=q10*par1*(1/(par2*exp((V-par3)/par4))+1/(par5+par6*exp(-(V-par7)/par8)))^(-1); % KHARCHE
dY(7,1) = (hs2_inf-hs2)/tau_hs2;

FsNa    = ((9.52e-02*exp(-6.3e-2*(V+34.4))/(1+1.66*exp(-0.225*(V+63.7))))+8.69e-2);
hs      = (1-FsNa)*hs1+FsNa*hs2;

% jhs not used
jhs_inf = hs2_inf;
tau_jhs = tau_hs2*1.5;
dY(41, 1) = (jhs_inf-jhs)/tau_jhs;

i_Nas = inas_block*g_Nas*ms^3*hs*Nao*(F1^2/(R2*T))*((exp((V-E_Na)*F1/(R2*T))-1)/(exp(V*F1/(R2*T))-1))*V; %pA/pF

%% INar

m  = Y(2);
h1 = Y(3);
h2 = Y(4);
jh = Y(40);

g_Na    = g_Nas*scaling_INa;
INa_act = -41;                     %mV
INa_In  = 65.1;                   %mV
sNa     = 4;

m_inf  = (1/(1+exp((V-INa_act)/(-sNa))));
tau_m   = q10*((0.6247/(0.832*exp(-0.335*(V+56.7))+0.627*exp(0.082*(V+65.01))))+0.04); %ms

tau_h1   = tau_hs1;
tau_h2   = tau_hs2;
dY(2, 1)= (m_inf-m)/tau_m;

h1_inf   = 1/(1+exp((V+INa_In)/4)); 
dY(3, 1) = (h1_inf-h1)/tau_h1;

h2_inf   = h1_inf;
dY(4, 1) = (h2_inf-h2)/tau_h2;
jh_inf = h2_inf;
tau_jh = tau_h2*1.5;
dY(40, 1) = (jh_inf-jh)/tau_jh;

FNa = ((9.52e-02*exp(-6.3e-2*(V+34.4))/(1+1.66*exp(-0.225*(V+63.7))))+8.69e-2);
h   = (1-FNa)*h1+FNa*h2;

i_Na = ina_block*g_Na*m^3*h*Nao*(F1^2/(R2*T))*((exp((V-E_Na)*F1/(R2*T))-1)/(exp(V*F1/(R2*T))-1))*V; %pA/pF

%% ICaL

d_L = Y(9);
f_L = Y(10);

g_Ca_L =0.07*scaling_ICaL; 
Cav1_2act= 10; 
slope1_2act = 4; 
Cav1_2In=36; %mV

alpha_d_L = (-0.01419*((V+35.0)/(exp(-(V+35)/2.5)-1))-0.04245*V/(exp(-0.208*V)-1)); %1/ms
beta_d_L  = (0.00571*(V-5)/(exp(0.4*(V-5))-1));  %1/ms
tau_d_L   = 1/(alpha_d_L+beta_d_L); %ms
d_L_inf   = 1/(1+exp(-(V+Cav1_2act)/slope1_2act));
dY(9,1)   = (d_L_inf-d_L)/tau_d_L;

alpha_f_L = (0.00312*(V+68)/(exp((V+68)/4)-1)); %1/ms
beta_f_L  = (0.025/(1+exp(-(V+68)/4))); %1/ms
tau_f_L   = 1/(alpha_f_L+beta_f_L); %ms
f_L_inf   = 1/(1+exp((V+Cav1_2In)/4.6));
dY(10,1)  = (f_L_inf-f_L)/tau_f_L;

i_Ca_L = ical_block*g_Ca_L*(f_L*d_L+(0.006/(1+exp(-(V+14.1)/6))))*(V-E_Ca); %pA/pF

%% ICaD

d_D = Y(11);
f_D = Y(12);

g_Ca_D=1.95*g_Ca_L; %nS/pF
%E_Ca_D=E_Ca_exp;
Cav1_3act = 22;
slope1_3act = 4.3;
Cav1_3In  = 48;

alpha_d_D = (-0.01419*((V+35.0)/(exp(-(V+35)/2.5)-1))-0.04245*V/(exp(-0.208*V)-1)); %1/ms
beta_d_D  = (0.00571*(V-5)/(exp(0.4*(V-5))-1));                                     %1/ms
tau_d_D   = 1/(alpha_d_D+beta_d_D);                                                 %ms
d_D_inf   = 1/(1+exp(-(V+Cav1_3act)/slope1_3act));
dY(11,1)  = (d_D_inf-d_D)/tau_d_D;

alpha_f_D = (0.00312*(V+68)/(exp((V+68)/4)-1)); %1/ms
beta_f_D  = (0.025/(1+exp(-(V+68)/4)));         %1/ms
tau_f_D   = 1/(alpha_f_D+beta_f_D);             %ms
f_D_inf   = 1/(1+exp((V+Cav1_3In)/5.4));
dY(12,1)  = (f_D_inf-f_D)/tau_f_D;

i_Ca_D = icaD_block*g_Ca_D*(f_D*d_D+(0.006/(1+exp(-(V+14.1)/6))))*(V-E_Ca); %pA/pF


%% ICaT
d_T = Y(13);
f_T = Y(14);

g_Ca_T   = icat_block*(0.0068/Cm*1000)*0.86; 

Cav3_1ac = 40;                    %mV 
Cav3_1In = 71;     %mV

alpha_d_T = 1.068*exp((V+26.3)/30);  %1/ms
beta_d_T  = 1.068*exp(-(V+26.3)/30); %1/ms
tau_d_T   = 1/(alpha_d_T+beta_d_T);  %ms
d_T_inf   = 1/(1+exp(-(V+Cav3_1ac)/6)); 
dY(13,1)  = (d_T_inf-d_T)/tau_d_T;

alpha_f_T = 0.0153*exp(-(V+71.7)/83.3); %1/ms
beta_f_T  = 0.015*exp((V+71.7)/15.38);  %1/ms
tau_f_T   = 1/(alpha_f_T+beta_f_T);     %ms
f_T_inf   = 1/(1+exp((V+Cav3_1In)/3.4));
dY(14,1)  = (f_T_inf-f_T)/tau_f_T;
% 
i_Ca_T = g_Ca_T*d_T*f_T*(V-E_Ca); %pA/pF

%%  IK1
gk1  = scaling_ik1*0.001525/Cm*1000; %nS/pF
i_K1 = gk1*(V-E_K)/(1+exp(0.07*(V+80.0))); %pA/pF


%% IKr
F_K_r = 0.0; 
P_af =  Y(15);
P_as =  Y(16);
P_i  =  Y(17);
P_a  =  (1-F_K_r)*P_af+F_K_r*P_as;

g_K_r = scaling_IKr*0.005/Cm*1000; %nS/pF

Pinac = -15;

P_af_inf = 1/(1+exp(-(V+23)/6.5));
tau_P_af = (0.84655/2.6./(0.0372.*exp((V-15)/20)+0.00096.*exp(-(V+40)/5)))+5; 
dY(15,1) = (P_af_inf-P_af)/tau_P_af;

P_as_inf = P_af_inf;
tau_P_as = (0.84655/(0.0042*exp(V/17)+0.00015*exp(-V/21.6)));  %ms
dY(16,1) = (P_as_inf-P_as)/tau_P_as;

tau_P_i  = 2;  %ms
P_i_inf  = 1/(1+exp((V+(-Pinac))/6.5));
dY(17,1) =(P_i_inf-P_i)/tau_P_i;

i_K_r = ikr_block*g_K_r*P_a*P_i*(V-E_K); %pA/pF

%% Ito and Isus

q = Y(18);
r = Y(19);

g_to  = (0.000491/Cm*1000)*scaling_Ito; %nS/pF
g_sus = scaling_Isus*0.0000665/Cm*1000;       %nS/pF

q_inf    = 1/(1+exp((V+59.37)/13.1));
tau_q    = (10.1+(65.17/((0.57*exp(-0.08*(V+49)))+0.000024*exp(0.1*(V+50.93)))));   %ms
dY(18,1) = (q_inf-q)/tau_q;

r_inf    = 1/(1+exp(-(V-10.93)/19.7));
tau_r    = (2.98+(15.59/(1.037*exp(0.09*(V+30.61))+0.369*exp(-0.12*(V+23.84))))); %ms
dY(19,1) = (r_inf-r)/tau_r;

i_to  = g_to*q*r*(V-E_K); %pA/pF
i_sus = g_sus*r*(V-E_K);  %pA/pF

%% Iks

xs = Y(23);
g_K_s  = scaling_IKs*0.000518/Cm*1000; %nS/pF

alpha_xs = (0.014/(1+exp(((-1)*(V-40))/9))); %1/ms
beta_xs  = (0.001*exp(((-1)*V)/45));         %1/ms
xs_inf   = (alpha_xs/(alpha_xs+beta_xs));
tau_xs   = (1/(alpha_xs+beta_xs));           %ms
dY(23,1) = (xs_inf-xs)/tau_xs;

i_K_s = g_K_s*xs^2*(V-E_K_s); %pA/pF

%% If

P = Y(20);
gh = 0.228*25;     % total conductance in Kharche
gh = (gh/Cm)*0.8; %nS/pF and 20% reduction to fit data

P_inf    = 1/(1+exp((V+97)/20));
tau_P    = 1.505./(exp(-0.0119*(V+590.3+14))+exp((V-55+14)/10)); %ms
dY(20,1) = (P_inf-P)/tau_P;

i_f_Na = 0.3833*gh*P*(V-E_K);  %pA/pF
i_f_K  = 0.3833*gh*P*(V-E_Na); %pA/pF
i_f = (i_f_Na+i_f_K)*if_block; %pA/pF

%% sodium_background_current

g_b_Na = 2*0.000058/Cm*1000; %nS/pF
i_b_Na = ibNa_block*g_b_Na*(V-E_Na); %pA/pF

%% potassium_background_current

g_b_K = 0.0000252/Cm*1000; %nS/pF
i_b_K = ibK_block*g_b_K*(V-E_K); %pA/pF

%% calcium_background_current

g_b_Ca = 0.5*0.0000132/Cm*1000; %nS/pF 
i_b_Ca = ibCa_block*g_b_Ca*(V-E_Ca); %pA/pF

%% sodium_potassium_pump

K_m_Na = 5.64;  %mM
K_m_K  = 0.621; %mM
i_p_max = 47.8*3; 

i_p = (i_p_max*(Nai/(K_m_Na+Nai))^3*(Ko/(K_m_K+Ko))^2*1.6/(1.5+exp(-(V+60)/40)))/Cm; %pA/pF

%% sodium calcium exchanger I NaCa from Paci model

KmCa   = 1.38;   % millimolar (in i_NaCa)
KmNai  = 87.5;   % millimolar (in i_NaCa)
Ksat   = 0.1;   % dimensionless (in i_NaCa)
gamma  = 0.35;  % dimensionless (in i_NaCa)

kNaCa		= 7000; % A_per_F (in i_NaCa)
kNaCa1 = kNaCa;   % A_per_F (in i_NaCa)
alpha		= 2.5;  % dimensionless (in i_NaCa)
i_NaCa = kNaCa1.*(exp(gamma.*V*F/(R*T))*Nai^3.0*Cao-exp((gamma-1.0).*V*F/(R*T))*Nao^3.0*Ca_sub*alpha)./((KmNai^3.0+Nao^3.0).*(KmCa+Cao)*(1.0+Ksat.*exp((gamma-1.0).*V*F/(R*T))));

%%  Ist
d_s = Y(21);
f_s = Y(22);
E_st = 10;
g_st = 0.017/Cm*1000*scaling_ist; %nS/pF 

alpha_d_s = 1/(0.15*exp(-V/11)+0.2*exp(-V/700)); %1/ms
beta_d_s  = 1/(16*exp(V/8)+15*exp(V/50));        %1/ms
tau_d_s   = 1/(alpha_d_s+beta_d_s);              %ms
d_s_inf   = alpha_d_s/(alpha_d_s+beta_d_s);
dY(21,1)  = (d_s_inf-d_s)/tau_d_s;

alpha_f_s = 1/(3100*exp(-V/13)+700*exp(-V/70));                    %1/ms
beta_f_s  = 1/(95*exp(-V/10)+50*exp(V/700))+2.29e-4/(1+exp(-V/5)); %1/ms
tau_f_s   = 1/(alpha_f_s+beta_f_s);                                %ms
f_s_inf   = alpha_f_s/(alpha_f_s+beta_f_s);
dY(22,1)  = (f_s_inf-f_s)/tau_f_s;

i_st = g_st*d_s*f_s*(V-E_st); %pA/pF

%% calcium_dependent_potassium_channel
g_SK  = 0.004; %nS/pF
m_cak = Y(24);

beta_k = 100;
n_SK   = 1.7;
Cac    = 0.000007; %mM
Car=((Ca_sub/Cac)^n_SK); %mM

m_inf_cak = Car/(1+Car);
tau_cak   = (10^(-3))/( beta_k*(1+Car)); %ms
dY(24,1)  = (m_inf_cak-m_cak)/tau_cak;

i_SK = g_SK*m_cak^2*(V-E_K); %pA/pF

%% background_muscarinic_potassium_channel_current
g_K_ACh=((0.0000)*(Ko^0.41));
i_K_ACh = (g_K_ACh*(Ki-(Ko*exp((((-1)*V)*F)/(R*T))))); %pA/pF

%% intracellular_calcium_dynamics

R1 = Y(25);
O1 = Y(26);
I1 = Y(27);
RI =(1-R1-O1-I1);

tau_dif_Ca = 0.4; %1/ms 
tau_tr     = 60;   %1/ms

ks   = 250;  %1/ms
P_up=0.04*pup_block;   %mM/ms
pumpkmf = 0.000246; %mM
pumpkmr = 3.29; %mM
nup = 2;
j_up      = P_up*((Cai/pumpkmf)^nup - (CaNSR/pumpkmr)^nup)/(1.0 + (Cai/pumpkmf)^nup + (CaNSR/pumpkmr)^nup);
j_tr      = ((CaNSR-CaJSR)/tau_tr); %mM/ms
j_SRCarel = jrel_block*ks*O1*(CaJSR-Ca_sub);   %mM/ms

MaxSR = 15;       %dimensionless
MinSR = 1;        %dimensionless
HSR   = 2.5;      %dimensionless
EC50_SR= 0.45;    %mM

koCa = 1.5;          %1/mM^2*ms
kom  = 0.06;         %1/ms
kiCa = 0.05;         %1/mM*ms
kim  = 0.005;        %1/mM*ms

kCaSR  = MaxSR-((MaxSR-MinSR)/(1+(EC50_SR/CaJSR)^HSR)); %dimensionless
koSRCa = (koCa/kCaSR); %1/mM^2*ms
kiSRCa = (kiCa*kCaSR); %1/mM*ms

dY(25,1) = (kim*RI-kiSRCa*Ca_sub*R1-(koSRCa*Ca_sub^2*R1-kom*O1)); %dR/dt
dY(26,1) = (koSRCa*Ca_sub^2*R1-kom*O1-(kiSRCa*Ca_sub*O1-kim*I1)); %dO/dt
dY(27,1) = (kiSRCa*Ca_sub*O1-kim*I1-(kom*I1-koSRCa*Ca_sub^2*RI)); %dI/dt

%% Ca+2 diffusion

j_Ca_dif =((Ca_sub-Cai)/tau_dif_Ca); %mM/ms



%% intracellular_ion_concentrations

%mM
CM_tot  = 0.045;
TC_tot  = 0.031;
TMC_tot = 0.062;
CQ_tot  = 10;

% 1/ms*mM
kf_TC  = 88.8;
kf_TMC = 227.7;
kf_TMM = 2.277;
kf_CM  = 227.7;
kf_CQ  = 0.534;

% 1/ms
kb_TC  = 0.446;
kb_TMC = 0.00751;
kb_TMM = 0.751;
kb_CM  = 0.542;
kb_CQ  = 0.445;

Mgi = 2.5; %mM

%mM
fTC  = Y(34);
fTMC = Y(35);
fTMM = Y(36);
fCMi = Y(37);
fCMs = Y(38);
fCQ  = Y(39);

%flux
delta_fTC  = (((kf_TC*Cai)*((1)-fTC))-(kb_TC*fTC));
delta_fTMC = (((kf_TMC*Cai)*((1)-(fTMC+fTMM)))-(kb_TMC*fTMC));
delta_fTMM = (((kf_TMM*Mgi)*((1)-(fTMC+fTMM)))-(kb_TMM*fTMM));
delta_fCMi = ((kf_CM*Cai)*((1)-fCMi))-(kb_CM*fCMi);
delta_fCMs = (((kf_CM*Ca_sub)*((1)-fCMs))-(kb_CM*fCMs));
delta_fCQ  = (((kf_CQ*CaNSR)*((1)-fCQ))-(kb_CQ*fCQ));

dY(28,1) = (-(i_f_Na+i_st+i_Nas+i_Na+i_b_Na+3*i_p+3*i_NaCa)*Cm)/(F*(V_i+V_sub)); %dNai/dt
dY(29,1) = (-(i_K_r+i_K_s+i_to+i_sus+i_SK+i_K1+i_f_K+i_K_ACh+(-2*i_p)+i_b_K-i_stim)*Cm/(F*(V_i+V_sub))); %dKi/dt
dY(30,1) = (j_Ca_dif*V_sub-j_up*V_NSR)/V_i -((CM_tot*delta_fCMi)+(TC_tot*delta_fTC)+(TMC_tot*delta_fTMC)); %dCa_i/dt
dY(31,1) = (-(i_Ca_L+i_Ca_T+i_b_Ca+i_Ca_D-2*i_NaCa)*Cm/(2*F)+((j_SRCarel)*V_JSR))/V_sub-(j_Ca_dif+(CM_tot*delta_fCMs)); %dCa_sub/dt
dY(32,1) = (j_up-((j_tr*V_JSR)/V_NSR));%dCaNSR/dt
dY(33,1) = (j_tr-(j_SRCarel+(CQ_tot*delta_fCQ))); %dCaJSR/dt

%% calcium_buffering

dY(34,1) = delta_fTC;  %dYfTC/dt
dY(35,1) = delta_fTMC; %dYfMC/dt
dY(36,1) = delta_fTMM; %dYfTMM/dt
dY(37,1) = delta_fCMi; %dYfCMi/dt
dY(38,1) = delta_fCMs; %dYfCMs/dt
dY(39,1) = delta_fCQ;  %dYfCQ/dt


%% Membrane potential
dY(1, 1) =  -(i_Na+i_Nas+i_Ca_L+i_Ca_D+i_Ca_T+i_K1+i_K_r+i_to+i_sus+ i_f+ i_K_s + i_b_Na+i_b_K+i_b_Ca+i_p+i_NaCa+ i_st+ i_SK+ i_K_ACh +number_of_macrophages*I_Gap_AV) +i_stim;

%% Output variables
INa   = i_Na;
INas  = i_Nas;
ICaL  = i_Ca_L;
ICaD  = i_Ca_D;
ICaT  = i_Ca_T;
Ik1   = i_K1;
Ikr   = i_K_r;
Ito   = i_to;
Isus  = i_sus;
If    = i_f;
IKs   = i_K_s;
IbNa  = i_b_Na;
IbK   = i_b_K;
IbCa  = i_b_Ca;
Ip    = i_p;
INaCa = i_NaCa;
Ist   = i_st;
Istim = i_stim;
ISK   = i_SK;
IKACh = i_K_ACh;
IGapPhi= I_Gap_Phi;
 IGapAV=I_Gap_AV;
IbPhi=Ib_Phi;
IK1Phi=IK1_Phi;
idiffPhi=i_diff_Phi;
IkurPhi=IKur_Phi;
IshkPhi=Ishk_Phi;

dati  = [INa, INas, ICaL , ICaD, ICaT, Ikr,Ik1, Ito, Isus, If, IbNa, IbK, IbCa, Ip, INaCa, Ist,...
    IKs, h1, h2,m, d_D, f_D, d_T,f_T, Istim, d_s, f_s, m_inf, tau_m,h1_inf, tau_h1, tau_h2,...
    d_D_inf, f_D_inf, tau_d_D,tau_f_D, d_T_inf, f_T_inf, tau_d_T,tau_f_T, ms, hs1, hs2,...
    ISK, m_cak , j_Ca_dif, j_up, j_tr, j_SRCarel, IKACh, Cai,Nai, Ki, Ca_sub, CaNSR,...
    CaJSR,O1, R1, I1, tau_P, xs, h, hs,r,q,P,IGapPhi, IGapAV,IbPhi,IK1Phi,idiffPhi,IkurPhi,IshkPhi];
end