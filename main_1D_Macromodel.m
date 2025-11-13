%%
close all
clear
clc



options = odeset('MaxStep',1, 'outputfcn', @odewbar);%,
tcouleurs=jet(2);
%************************************************
% Initial conditions for the basal cell model
%**********************************************
Nbrecell=24;
%load('Init_State.mat');
Y=zeros(50,Nbrecell);
%Yc=zeros(41,Nbrecell);

Y(1,1:Nbrecell)= -60.1; % potential
Y(2,1:Nbrecell)=0.135; % m activation gate INa
Y(3,1:Nbrecell)=0.03123; %h1 inactivation gate INa
Y(4,1:Nbrecell)=1;%4.894; %h2 inactivation gate INa
Y(5,1:Nbrecell)=0.135; % ms activation gate INas
Y(6,1:Nbrecell)=0.03123; %h1s inactivation gate INas
Y(7,1:Nbrecell)=1;%4.894; %h2s inactivation gate INas

Y(8,1:Nbrecell)  = 0.02248; % d_L activation gate ICaL Cav1_2
Y(9,1:Nbrecell)=0.43; %f_L inactivation gate ICaL Cav1_2

Y(10,1:Nbrecell)=0.02248; %d_D activation gate ICaL Cav1_3
Y(11,1:Nbrecell)=0.531; %f_D inactivation gate ICaL Cav1_3

Y(12,1:Nbrecell)=0.02217; %d_T activation gate ICaT
Y(13,1:Nbrecell)=0.06274; % f_T inactivation gate ICaT

Y(14,1:Nbrecell)=0.002359413; %P_af activation fast gate Ikr
Y(15,1:Nbrecell)=0.09102082; %P_as activation slow gate Ikr
Y(16,1:Nbrecell)=0.9977152; % P_i inactivation gate Ikr

Y(17,1:Nbrecell)=0.06609; % q inactivation gate Ito
Y(18,1:Nbrecell)=0.05733; % r activation gate Ito and I_sus

Y(19,1:Nbrecell)=0.007645; %P hyperpolarisation_activated_current_y_gate

Y(20,1:Nbrecell)=0.0015225; %d_s sustained_inward_current_d_gate
Y(21,1:Nbrecell)=0.283; %f_s sustained_inward_current_f_gate

Y(22,1:Nbrecell)=0; %xs gate slow_delayed_rectifying_potassium_current

Y(23,1:Nbrecell)=0.99979196;%0; %m_cak

Y(24,1:Nbrecell)=0.64484881; %R1
Y(25,1:Nbrecell)=0.000012240116; %O1
Y(26,1:Nbrecell)=0.0000067409369; %I1

Y(27,1:Nbrecell)=7; % Nai Update to the SS value
Y(28,1:Nbrecell)=159.43881; %Ki
Y(29,1:Nbrecell)=0.0000326892; %Cai
Y(30,1:Nbrecell)=0.0000534741; %Ca_sub
Y(31,1:Nbrecell)=4.5887879; %CaNSR
Y(32,1:Nbrecell)=3.431262; %CaJSR

Y(33,1:Nbrecell)=0.45367981;  %fTC
Y(34,1:Nbrecell)=0.93724623;  %fMC
Y(35,1:Nbrecell)=0.05542706;  %fTMM
Y(36,1:Nbrecell)=0.6364781;   %fCMi
Y(37,1:Nbrecell)=0.29785784;  %fCMs
Y(38,1:Nbrecell)=0.84628533;  %fCQ
Y(39,1:Nbrecell)=1;            %jh
Y(40,1:Nbrecell)=1;            %jhs
Y(41,1:Nbrecell)=-50.6634; %potential Macrophage mV
Y(42,1:Nbrecell)=1;% CO
Y(43,1:Nbrecell)=0; % C1
Y(44,1:Nbrecell)=0;% C2
Y(45,1:Nbrecell)=0;% C3
Y(46,1:Nbrecell)=0;% C4
Y(47,1:Nbrecell)=0; % O
Y(48,1:Nbrecell)=0;% I
Y(49,1:Nbrecell)=0;% iKur.R
Y(50,1:Nbrecell)=0;% S


Y=Y(:)';

%**************************************************************************

%% Integration method
model = @AVNmodel_Macromodel_1D;




%****************************************************
% fibre Matrice
%************************************************
%save("Init_1DCtrl_from40bpm550coupling12iter.mat","Y");
%load("Init_1DCtrl.mat");
 % Y1=Y(1:42);
 % Y=repmat(Y1,1,Nbrecell);
%  2700 2800
facteur_echelle_couplage=3000;
facteur_echelle_couplage_Macro=7;
nsec =20; %s
CL = 1000*nsec;
G_gapm=facteur_echelle_couplage_Macro*1;%40;   % macrophage conductance nS
G_gap=facteur_echelle_couplage*40;     % gap junction coupling nS
Cm_Phi=18.32;%27.9;
inter_macro=4;
group_macro=[0,6,7,8];
Macro_Type=12;
%***********************************************************
bpm_stim =350; %      stimulation  frequency BPM
flag_stim=1;   % flag for activated stimulation 
   i_stim_Amplitude 		=0.2*1e3;   %nanoampere (in stim_mode)
        %i_stim_Amplitude_Macro=0.6*1e3;
        i_stim_End 				= CL;   % msecond (in stim_mode)
        i_stim_PulseDuration	=5;   % msecond (in stim_mode)
        i_stim_Start 			= 0;   % msecond (in stim_mode)
        i_stim_frequency 		=bpm_stim;   % beat per_minut (in stim_mode)
        i_stim_Period 			=60*1e3/i_stim_frequency; %ms
        NbreCellules_a_stimule=2;

%
Nbre_episodes=1;
Temps_Application_Macro=10000;
dateSimulation='06112025';
ActivShift=0;
i_diff_Phi=0;

%*****************************************************************************

 chemin=['F:\ResultatsSimulationMacrophages\Simulation_1D_',dateSimulation,'\Resultat_Ncells_',num2str(Nbrecell),'_typeM_',num2str(Macro_Type),'_gap_',...
   num2str(G_gap),'_Macro_',num2str(G_gapm),'_Inter_Macro_',num2str(inter_macro),'_ActivShift',num2str(ActivShift),'_GroupeMacro_',num2str(group_macro),...
   '_Stim_',num2str(flag_stim),'_BPM_',num2str(bpm_stim),...
   '_NbrecellulesStim_',num2str(NbreCellules_a_stimule),...
    '_i_stim_PulseDuration_',num2str(i_stim_PulseDuration),...
   ' _i_stim_Amplitude_',num2str( i_stim_Amplitude) ,'_time_',num2str(CL),'_sec\'];  % num2str(G_gapm)

  test_dir=exist(chemin,'dir');
    if test_dir==0
        mkdir(chemin)
    end

%**************************************************************
% Cell dimension
%******************************************************************
Cm      =  22;     %pF
Rcell   = 0.00005; %dm
Cm_surf = 1;       %uF/cm^2
Lcell = (1/(2*pi*Rcell))*(Cm/(10^8*Cm_surf)-(2*pi*(Rcell)^2)); %dm % 65 um
deltax=65; % deta distance um
L=deltax*Nbrecell;
%distance=0:dx;L;
%%
for episode=1:Nbre_episodes % run by episode of CL  time length msec
   
    tic
    Const=[bpm_stim,flag_stim,G_gap,G_gapm,Cm_Phi,CL,deltax,...
        i_stim_End, i_stim_PulseDuration,i_stim_Start,i_stim_Amplitude,...
        i_stim_frequency,i_stim_Period, NbreCellules_a_stimule,inter_macro,...
        Temps_Application_Macro,Macro_Type,group_macro,ActivShift,Nbrecell,i_diff_Phi];

    [t ,Yout] = ode15s(@(x,y) model(x, y,Const),[0 CL],Y , options);
%%
    IgapAV  =zeros(Nbrecell,size(Yout,1));
    IgapMyo_Macro=zeros(Nbrecell,size(Yout,1));
    for i= 1:size(Yout,1)
        [~, dati]    = model(t(i), Yout(i,:),Const);
        IgapAV(:,i)=dati(1:Nbrecell);
        IgapMyo_Macro(:,i)=dati(25:Nbrecell*2);

    end
%
%****************************************************
% *******
%
%**********************************************************************


Vm=zeros(size(Yout,1),Nbrecell);
k=1;
for i=1:50:size(Yout,2)
Vm(:,k)=Yout(:,i);
k=k+1;
end
           
Y=Yout(end,:);

elapsedTime = toc;
disp(elapsedTime/60)
end
%%
figure(episode)


imagesc(t,(1:Nbrecell).*deltax,Vm')
xlabel('time msec');ylabel('distance um')

%******************************************
figure(episode+Nbrecell)
for cell=1:Nbrecell-18
subplot(6,1,cell)
plot(t,Vm(:,cell))
yyaxis  right
plot(t,IgapAV(cell,:))
hold on
plot(t,IgapMyo_Macro(cell,:),'g')
%plot([t(1) t(end)],[0 0])
yyaxis left
ylabel(['Cell ',num2str(cell)])

end
%***********************************
figure(episode+Nbrecell*2)
for cell=7:Nbrecell-12
subplot(6,1,cell-6)
plot(t,Vm(:,cell))
yyaxis  right
plot(t,IgapAV(cell,:))
hold on
plot(t,IgapMyo_Macro(cell,:),'g')
%plot([t(1) t(end)],[0 0])
yyaxis left
ylabel(['Cell ',num2str(cell)])

end
%**********************************
figure(episode+Nbrecell*3)
for cell=13:Nbrecell-6
subplot(6,1,cell-12)
plot(t,Vm(:,cell))
yyaxis  right
plot(t,IgapAV(cell,:))
hold on
plot(t,IgapMyo_Macro(cell,:),'g')
%plot([t(1) t(end)],[0 0])
yyaxis left
ylabel(['Cell ',num2str(cell)])
end
%***************************************
figure(episode+Nbrecell*4)
for cell=19:Nbrecell
subplot(6,1,cell-18)
plot(t,Vm(:,cell))
yyaxis  right
plot(t,IgapAV(cell,:))
hold on
plot(t,IgapMyo_Macro(cell,:),'g')
%plot([t(1) t(end)],[0 0])
yyaxis left
ylabel(['Cell ',num2str(cell)])
end
%****************************************
figure(episode+Nbrecell*5)
imagesc(t,(1:Nbrecell),IgapAV)

%****************************************
figure(episode+Nbrecell*6)
imagesc(t,(1:Nbrecell),IgapMyo_Macro)
%*********************************************************************
%
%***********************************************************************

%***************************************************************************     
dx=0.0065028;
%Nbrecell=Ncells;
indextime=find(t<10000,1,"last");
BpmTable=zeros(Nbrecell,4);
for cell=1:Nbrecell
 Tf=islocalmax(Vm(:,cell),'minProminence',45,'MinSeparation',100,'SamplePoints',t);
Tmdp=islocalmin(Vm(:,cell),'minProminence',45,'MinSeparation',150,'SamplePoints',t);
Events{cell}=Tf;
EventsMDP{cell}=Tmdp;




Nbre1=sum(Tf(1:indextime));
Nbre2=sum(Tf(indextime:end));
BpmTable(cell,1)=Nbre1;
BpmTable(cell,2)=Nbre2;
BpmTable(cell,3)=60*Nbre1/10;
BpmTable(cell,4)=60*Nbre2/10;


end
%********************************************************************************
% recherche du seuil et 90% repolarisation
%**********************************************************************************
clc
dt=mean(diff(t));
dvdt=Vm;
for cell=1:Nbrecell
    
dvdt(:,cell)=[diff(Vm(:,cell));0];

 Tf=islocalmax(Vm(:,cell),'minProminence',45,'MinSeparation',100,'SamplePoints',t);
Tmdp=islocalmin(Vm(:,cell),'minProminence',45,'MinSeparation',150,'SamplePoints',t);
%Tdvdt=islocalmax(dvdt(:,cell),'minProminence',0.5,'MinSeparation',150,'SamplePoints',t);
Events{cell}=Tf;
EventsMDP{cell}=Tmdp;
   
Tmax=Events{cell};
Tmdp=EventsMDP{cell};

indexTmax=find(Tmax);
indexTmdp=find(Tmdp);

Fenetre_detection=30;
TSeuil=false(size(Vm,1),1);
T90=false(size(Vm,1),1);
for i=1:length(indexTmax)-1
Amp=Vm(indexTmax(i),cell);
  [s, ~]=Recher_Start_Stop(dvdt(:,cell),indexTmax(i),0.015,Fenetre_detection);
  if s==0
      s=1;
  end


 TSeuil(s,1)=true;
 Inmdp=find(indexTmdp>indexTmax(i),1,'first');
 
 Amp90=Amp -(Amp-Vm(indexTmdp(Inmdp),cell))*0.90;
index90=indexTmax(i)+find(Vm(indexTmax(i):indexTmdp(Inmdp),cell)>Amp90,1,'last');
T90(index90,1)=true;

end

T90repol{cell}=T90;
SeuilPa{cell}=TSeuil;



end
%***************************************************************************
figure(episode+Nbrecell*7)
imagesc(Vm')
hold on
for cell=1:Nbrecell
Tpa=Events{cell};
Tmdp=EventsMDP{cell};
T90=T90repol{cell};
TSeuil=SeuilPa{cell};
%*******************************
indextimepa=find(Tpa);
indextimeMDP=find(Tmdp);
indextimeT90=find(T90);
indextimeSeuil=find(TSeuil);
%**********************************
pamax=indextimepa;
pamax(:)=cell;
pamdp=indextimeMDP;
pamdp(:)=cell;
paseuil=indextimeSeuil;
paseuil(:)=cell;
paT90=indextimeT90;
paT90(:)=cell;

plot(indextimepa,pamax,'*r','MarkerSize',10)
plot(indextimeT90,paT90,'og','MarkerSize',10)
plot(indextimeSeuil,paseuil,'+r','MarkerSize',10)
plot(indextimeMDP,pamdp,'>c','MarkerSize',10)
end
hold off
%****************************************************************

for i=1:24
    nc=i;

    figure(nc)


    Tpa=Events{nc};
    Tmdp=EventsMDP{nc};
    T90=T90repol{nc};
    TSeuil=SeuilPa{nc};
    %*******************************
    indextimepa=find(Tpa);
    indextimeMDP=find(Tmdp);
    indextimeT90=find(T90);
    indextimeSeuil=find(TSeuil);
   
    plot(t,Vm(:,nc))
    hold on
    plot(t(indextimepa),Vm(indextimepa,nc),'*r','MarkerSize',10)
    plot(t(indextimeT90),Vm(indextimeT90,nc),'og','MarkerSize',10)
    plot(t(indextimeSeuil),Vm(indextimeSeuil,nc),'+r','MarkerSize',10)
    plot(t(indextimeMDP),Vm(indextimeMDP,nc),'>c','MarkerSize',10)
    hold off
   Ioff=find(t(indextimeSeuil)<10000);
   Ion=find(t(indextimeSeuil)>=10000);
 APD90off{nc}=t(indextimeT90(Ioff))-t(indextimeSeuil(Ioff));
 APD90on{nc}=t(indextimeT90(Ion))-t(indextimeSeuil(Ion));
 
tableduree=[t(indextimeSeuil),t(indextimeT90),t(indextimeT90)-t(indextimeSeuil)];
end
%%
 Pool_APD90off=[];Pool_APD90on=[];
for i=1:24

tmp_APD= APD90off{i} ;
Pool_APD90off=[Pool_APD90off; tmp_APD];

tmp_APD= APD90on{i} ;
Pool_APD90on=[Pool_APD90on ;tmp_APD];
 
end


%


%********************************************************************************
Nc1=5;Nc2=10;
Tf=Events{Nc1};
Tl=Events{Nc2};
indexT1=find(Tf(1:indextime));
indexT2=find(Tl(1:indextime));


vconduction1=zeros(length(indexT1)-1,1);

distance=Nc2*dx-Nc1*dx;

Tf_pamax=t(indexT1).*1e-3;
Tl_pamax=t(indexT2).*1e-3;

for i=1:length(Tf_pamax)
Tpaf=Tf_pamax(i);
indexl=find(Tl_pamax>Tpaf,1,"first");
Tpal=Tl_pamax(indexl);
duree=Tpal-Tpaf;
vconduction1(i)=distance/duree;
duree1(i)=duree;
end

indexT1=find(Tf(indextime:end));
indexT2=find(Tl(indextime:end));
vconduction2=zeros(length(indexT1)-1,1);

Tf_pamax=t(indexT1).*1e-3;
Tl_pamax=t(indexT2).*1e-3;
for i=1:length(Tf_pamax)
Tpaf=Tf_pamax(i);
indexl=find(Tl_pamax>Tpaf,1,"first");
Tpal=Tl_pamax(indexl);
duree=Tpal-Tpaf;
vconduction2(i)=distance/duree;
duree2(i)=duree;
end
%******************************************************************************************
% Affiche Resultat
%******************************************************************************************
figure(episode+Nbrecell*8)
subplot(1,2,1)
histogram(Pool_APD90off,30) 
hold on
histogram(Pool_APD90on,30) 
subplot(1,2,2)
histogram(vconduction1,30) 
hold on
histogram(vconduction2,30) 
%*******************************************************************************************
% sauver resultat analyse
%******************************************************************************************
save([chemin,'simul',num2str(Nbrecell),'cells_couple',num2str(G_gap),'_',...
    num2str(episode),'_',num2str(bpm_stim),'_couplemacro0_stim_',...
    num2str(flag_stim),'.mat'],'-mat');