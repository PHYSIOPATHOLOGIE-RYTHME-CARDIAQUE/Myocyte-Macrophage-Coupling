close all
clear
clc
perStim=[0,450,550];
for i=1:1
    flag_perstim=1;

    mainAVNMacro(perStim(i))

end
function mainAVNMacro(pst)

%***************************************************************************

    linecolors = jet(4);
    model = @AVNMacromodel_new;
    nsec = 20; %s
    CL= 1000*nsec;
    %flagsim ='test';
    session=150;
    resultat_sauv_affiche=0;
    gapstep=0.05;
    G_gap=0;
    i_ke=1;
    ResultVm=zeros(1,10);
    %K_e_Phi=5.4;
    Cm_Phi=18.32;%27.9;
    %  tableCouleur=jet(session);
    %% Integration method
    number_of_macrophages=1;
    liste_K={5.4,10,50,120};
    K_e_Phi=liste_K{1};

    liste_Sym={'o','s','*','<'};
    bpm_stim =pst;%fps*round(174+174*pst);%212; %CL (ms)

    Choix_Type_Phi=0;

    i_stim_Amplitude 		=0.6*1e3;   %nanoampere (in stim_mode)
    i_stim_Amplitude_Macro=0.6*1e3;
    i_stim_End 				= CL*1000;   % msecond (in stim_mode)
    i_stim_PulseDuration	=1;   % msecond (in stim_mode)
    i_stim_Start 			= 0;   % msecond (in stim_mode)
    i_stim_frequency 		=bpm_stim;   % beat per_minut (in stim_mode)
    stim_flag 	=0;
    stim_flagM=0;
    % if sim==200
    % stim_flag 				=0;   % dimensionless (in stim_mode)
    % else
    %    stim_flag 				=1;   % dimensionless (in stim_mode)
    % end
    % stim_flagM= 0;
    i_stim_Period 			=60*1e3/i_stim_frequency; %ms
    minPerioddetect=0.1;
    minAmpdetect=5;
    %Rate=zeros(1,2);
    %par_block=1;
    tic
    options = odeset('MaxStep',1, 'outputfcn', @odewbar);

    %path=pwd();
    %dir_result=[ path,'\ResultatsDXround_typeMPhi',num2str( Choix_Type_Phi),'StimMacro',num2str( stim_flagM),...
    %   'StimMyo',num2str( stim_flag),'gapstep005nSsession',num2str(session),'capaM',num2str(Cm_Phi),'\'];
    dir_result=[ 'F:\ResultatsSimulationMacrophages\Resultats04112025_typeMPhi',num2str( Choix_Type_Phi),...
        'GgapStep',num2str(gapstep),'StimMacro',num2str( stim_flagM),...
        'StimMyo',num2str( stim_flag),'perStim',num2str( bpm_stim),'AmpSt',num2str(i_stim_Amplitude),...
        'session',num2str(session),'capaM',num2str(Cm_Phi),'\'];
    test_dir=exist(dir_result,'dir');
    %%
    if test_dir==0
        mkdir(dir_result)
    end


    %tcouleurs=jet(2);
    %% Initial conditions for the basal cell model
   

    % Y(1) = -60.1; % potential
    % Y(2) = 0.135; % m activation gate INa
    % Y(3) = 0.03123; %h1 inactivation gate INa
    % Y(4) = 1;%4.894; %h2 inactivation gate INa
    % Y(5) = 0.135; % ms activation gate INas
    % Y(6) = 0.03123; %h1s inactivation gate INas
    % Y(7) = 1;%4.894; %h2s inactivation gate INas
    %
    % Y(9)  = 0.02248; % d_L activation gate ICaL Cav1_2
    % Y(10) = 0.43; %f_L inactivation gate ICaL Cav1_2
    %
    % Y(11) = 0.02248; %d_D activation gate ICaL Cav1_3
    % Y(12) = 0.531; %f_D inactivation gate ICaL Cav1_3
    %
    % Y(13) = 0.02217; %d_T activation gate ICaT
    % Y(14) = 0.06274; % f_T inactivation gate ICaT
    %
    % Y(15) = 0.002359413; %P_af activation fast gate Ikr
    % Y(16) = 0.09102082; %P_as activation slow gate Ikr
    % Y(17) = 0.9977152; % P_i inactivation gate Ikr
    %
    % Y(18) = 0.06609; % q inactivation gate Ito
    % Y(19) = 0.05733; % r activation gate Ito and I_sus
    %
    % Y(20) = 0.007645; %P hyperpolarisation_activated_current_y_gate
    %
    % Y(21) = 0.0015225; %d_s sustained_inward_current_d_gate
    % Y(22) = 0.283; %f_s sustained_inward_current_f_gate
    %
    % Y(23) = 0; %xs gate slow_delayed_rectifying_potassium_current
    %
    % Y(24) = 0.99979196;%0; %m_cak
    %
    % Y(25) = 0.64484881; %R1
    % Y(26) = 0.000012240116; %O1
    % Y(27) = 0.0000067409369; %I1
    %
    % Y(28) = 7; % Update to the SS value
    % Y(29) = 159.43881; %Ki
    % Y(30) = 0.0000326892; %Cai
    % Y(31) = 0.0000534741; %Ca_sub
    % Y(32) = 4.5887879; %CaNSR
    % Y(33) = 3.431262; %CaJSR
    %
    % Y(34) = 0.45367981;  %fTC
    % Y(35) = 0.93724623;  %fMC
    % Y(36) = 0.05542706;  %fTMM
    % Y(37) = 0.6364781;   %fCMi
    % Y(38) = 0.29785784;  %fCMs
    % Y(39) = 0.84628533;  %fCQ
    % Y(40) = 1;            %jh
    % Y(41) = 1;            %jhs
    % Y(42)=1;
    % Y(43)=0;
    % Y(44)=0;
    % Y(45)=0;
    % Y(46)=0;
    % Y(47)=1;
    % Y(48)=0;
    % Y(49)=0;
    % Y(50)=0;
    % Y(51)=-51; %V membrane Macrophages mV



    for sim=1:session


       
        Init=load('initMacro.mat'); % initial conditions after 1 minutes run.
        Y=Init.Y;

        Const=[G_gap,Choix_Type_Phi,K_e_Phi,i_stim_Amplitude,i_stim_End,i_stim_PulseDuration,i_stim_Start, i_stim_frequency,...
            stim_flag, stim_flagM,i_stim_Period,number_of_macrophages, i_stim_Amplitude_Macro,Cm_Phi];





        %flagsim2 = ['lastMDP_stim', num2str(bpm_stim),'bpm_'];
        %flagsim  = ['data_stim', num2str(bpm_stim),'bpm_'];
        % flagsim2 = ['lastEnd_ICaDBlock_'];
        % flagsim  = ['dataEndICaDBlock_'];
        %flag_save = 0;
        % load CI_old.mat
        % Y=Yfinal;
        % % loading of steady state IC
        %lload CIss_lastMDP_300s
        %Y=IC';
        % Y(1)=Y(1);
        [t,Yc] = ode15s(@(x,y) model(x, y, Const),[0 CL], Y , options);
        Vm   = Yc(:,1);
        dVm  = [0; diff(Vm)./diff(t)];
        %Yfinal = Yc(end,:);
        VmMphi=Yc(:,51);
        dVmPhi=[0;diff(VmMphi)./diff(t)];
        %***************************************************************************************************
        INa   =zeros(1,size(Yc,1));INas  =zeros(1,size(Yc,1)); ICaL  = zeros(1,size(Yc,1)); ICaD  = zeros(1,size(Yc,1));
        ICaT  = zeros(1,size(Yc,1)); Ikr   = zeros(1,size(Yc,1));Ik1   = zeros(1,size(Yc,1));Ito   = zeros(1,size(Yc,1));Isus  = zeros(1,size(Yc,1));
        If    = zeros(1,size(Yc,1)); IbNa  = zeros(1,size(Yc,1)); IbK   = zeros(1,size(Yc,1));IbCa  = zeros(1,size(Yc,1));
        Ip    = zeros(1,size(Yc,1));INaCa = zeros(1,size(Yc,1));Ist   = zeros(1,size(Yc,1));IKs   = zeros(1,size(Yc,1));h1    =zeros(1,size(Yc,1));
        h2    =zeros(1,size(Yc,1));    m     =zeros(1,size(Yc,1));    d_D   =zeros(1,size(Yc,1));    f_D   =zeros(1,size(Yc,1));
        d_T   =zeros(1,size(Yc,1));    f_T   =zeros(1,size(Yc,1));    Istim =zeros(1,size(Yc,1));    d_s   =zeros(1,size(Yc,1));
        f_s   =zeros(1,size(Yc,1));    m_inf  = zeros(1,size(Yc,1));  tau_m  = zeros(1,size(Yc,1));  h1_inf = zeros(1,size(Yc,1));
        tau_h1 = zeros(1,size(Yc,1));  tau_h2 = zeros(1,size(Yc,1));  d_D_inf  =zeros(1,size(Yc,1));
        f_D_inf  =zeros(1,size(Yc,1));  tau_d_D  =zeros(1,size(Yc,1));  tau_f_D  =zeros(1,size(Yc,1));  d_T_inf  =zeros(1,size(Yc,1));
        f_T_inf  =zeros(1,size(Yc,1));tau_d_T  =zeros(1,size(Yc,1)); tau_f_T  =zeros(1,size(Yc,1));
        ms   =zeros(1,size(Yc,1));hs1  =zeros(1,size(Yc,1)); hs2  =zeros(1,size(Yc,1)); ISK    = zeros(1,size(Yc,1)); m_cak  = zeros(1,size(Yc,1));
        j_Ca_dif  = zeros(1,size(Yc,1)); j_up      = zeros(1,size(Yc,1)); j_tr      = zeros(1,size(Yc,1));
        j_SRCarel = zeros(1,size(Yc,1));IKACh  = zeros(1,size(Yc,1)); Cai = zeros(1,size(Yc,1));  Nai = zeros(1,size(Yc,1)); Ki  = zeros(1,size(Yc,1));
        Ca_sub = zeros(1,size(Yc,1));CaNSR  = zeros(1,size(Yc,1));CaJSR  = zeros(1,size(Yc,1));
        h  = zeros(1,size(Yc,1)); hs  = zeros(1,size(Yc,1));r=zeros(1,size(Yc,1)); q=zeros(1,size(Yc,1)); P=zeros(1,size(Yc,1));
        I_Gap_Phi=zeros(1,size(Yc,1));  I_Gap_AV=zeros(1,size(Yc,1)); Ib_Phi=zeros(1,size(Yc,1)); IK1_Phi=zeros(1,size(Yc,1));i_diff_Phi=zeros(1,size(Yc,1));
        IKur_Phi=zeros(1,size(Yc,1)); Ishk_Phi=zeros(1,size(Yc,1));
        %***************************************************************************************************
        %%
        for i= 1:size(Yc,1)
            [~, dati]    = model(t(i), Yc(i,:),Const);
            INa(i)   = dati(1);
            INas(i)  = dati(2);
            ICaL(i)  = dati(3);
            ICaD(i)  = dati(4);
            ICaT(i)  = dati(5);
            Ikr(i)   = dati(6);
            Ik1(i)   = dati(7);
            Ito(i)   = dati(8);
            Isus(i)  = dati(9);
            If(i)    = dati(10);
            IbNa(i)  = dati(11);
            IbK(i)   = dati(12);
            IbCa(i)  = dati(13);
            Ip(i)    = dati(14);
            INaCa(i) = dati(15);
            Ist(i)   = dati(16);
            IKs(i)   = dati(17);
            h1(i)    = dati (18);
            h2(i)    = dati (19);
            m(i)     = dati (20);
            d_D(i)   = dati (21);
            f_D(i)   = dati (22);
            d_T(i)   = dati (23);
            f_T(i)   = dati (24);
            Istim(i) = dati (25);
            d_s(i)   = dati (26);
            f_s(i)   = dati (27);

            %gates INas
            m_inf(i)  = dati(28);
            tau_m(i)  = dati(29);
            h1_inf(i) = dati(30);
            tau_h1(i) = dati(31);
            tau_h2(i) = dati(32);

            %Gates of ICaD
            d_D_inf (i) = dati (33);
            f_D_inf (i) = dati (34);
            tau_d_D (i) = dati (35);
            tau_f_D (i) = dati (36);

            %Gates of ICaT
            d_T_inf (i) = dati (37);
            f_T_inf (i) = dati (38);
            tau_d_T (i) = dati (39);
            tau_f_T (i) = dati (40);

            %INas
            ms (i)  = dati (41);
            hs1 (i) = dati (42);
            hs2 (i) = dati (43);

            %ISK
            ISK (i)   = dati(44);
            m_cak (i) = dati(45);

            %Ca flux
            j_Ca_dif (i) = dati(46);
            j_up (i)     = dati(47);
            j_tr (i)     = dati(48);
            j_SRCarel (i)= dati(49);

            %IKACh
            IKACh (i) = dati(50);

            %intracellular currents
            Cai(i) = dati(51);
            Nai(i) = dati(52);
            Ki(i)  = dati(53);

            %calcium currents
            Ca_sub(i) = dati(54);
            CaNSR(i)  = dati(55);
            CaJSR(i)  = dati(56);
            h(i)  = dati(62);
            hs(i)  = dati(63);

            % Ito/IKur currents
            r(i)=dati(64);
            q(i)=dati(65);

            %If current
            P(i)=dati(66);
            I_Gap_Phi(i)=dati(67);
            I_Gap_AV(i)=dati(68);
            Ib_Phi(i)=dati(69);
            IK1_Phi(i)=dati(70);
            i_diff_Phi(i)=dati(71);
            IKur_Phi(i)=dati(72);
            Ishk_Phi(i)=dati(73);


        end
        t = t/1e3; % ms to s

        % result       = [INa;INas; ICaL; ICaD; ICaT; Ikr; Ik1; Ito ; Isus; If; IbNa; IbK; IbCa; Ip;...
        %     INaCa; Ist; IKs; h1; h2;m; d_D; f_D; d_T;f_T; Istim; d_s; f_s; m_inf; tau_m ; h1_inf;...
        %     tau_h1; tau_h2; d_D_inf; f_D_inf; tau_d_D; tau_f_D; d_T_inf; f_T_inf; tau_d_T;...
        %     tau_f_T; ms; hs1; hs2; ISK; m_cak; j_Ca_dif; j_up; j_tr; j_SRCarel; IKACh; r;q; P];
         I_tot_myo = sum([INa;INas; ICaL; ICaD; ICaT; Ikr; Ik1; Ito; Isus; If; IbNa; IbK; ...
             IbCa; Ip; INaCa; Ist; IKs; Istim; ISK; IKACh;I_Gap_AV],1);

        I_tot_macro=sum([Ib_Phi;IK1_Phi;i_diff_Phi;IKur_Phi;Ishk_Phi;I_Gap_Phi],1);

        mat_correnti = [t';Vm';VmMphi';INa;INas; ICaL; ICaD; ICaT; Ikr; Ik1; Ito; Isus; If; IbNa; IbK; ...
            IbCa; Ip; INaCa; Ist; IKs; Istim; ISK; IKACh;  I_Gap_Phi;...
            I_Gap_AV;Ib_Phi;IK1_Phi;i_diff_Phi;IKur_Phi;Ishk_Phi;dVm';dVmPhi'; I_tot_myo;I_tot_macro];




        PAR_Results=array2table(mat_correnti', 'VariableNames',{'time','Vm','VmMphi','INa','INas',...
            'ICaL',' ICaD','ICaT','Ikr','Ik1','Ito','Isus','If','IbNa','IbK','IbCa','Ip','INaCa','Ist','IKs',...
            'Istim','ISK','IKACh','IGap_Phi','IGap_AV','Ib_Phi','IK1_Phi','i_diff_Phi','IKur_Phi','Ishk_Phi',  'dVm','dVmMphi','inetmyo','inetmacro'});
        NomRes=[dir_result,'MacrophageMyo_Resultats','_type',num2str(Choix_Type_Phi),'session',num2str(sim),'.csv'];
        writetable(  PAR_Results,NomRes,'Delimiter','tab')


        %data= [string(datetime('today')),'_macro'];
       % file = [data,flagsim,num2str(nsec),'s.mat'];
       % file_IC = [data,flagsim2,num2str(nsec),'s.mat'];
        %[mdp index_mdp] = findpeaks(-Vm);
        toc
        % % %
        % if flag_save == 1
        %     %     Yfinal = Yc(index_mdp(end),:);
        %     IC = Yfinal;
        %     save(fullfile(file_IC), 'IC')
        %     save(fullfile(file),'t','Vm', 'Cai', 'Ca_sub', 'CaNSR', 'CaJSR', 'Nai', 'I_tot', 'result' ,'Yfinal')
        % end

        %% Plot








        TF = islocalmax(Vm,'MinProminence', minAmpdetect,'MinSeparation', minPerioddetect,'SamplePoints',t);
        TB= islocalmin(Vm,'MinProminence', minAmpdetect,'MinSeparation', minPerioddetect,'SamplePoints',t);
        % TFd=islocalmax(dVm,'MinProminence',4,'MinSeparation', minPerioddetect,'SamplePoints',t);
        %  indexmaxPA=find(TF);
        %  indexMDPPA=find( TB);


        %dVm
        %derivVm  = [diff(Vm);0];
        %%
        %d2Vm  = [0; diff(dVm)./diff(t)];


        figure('units','normalized','outerposition',[0 0 1 1])



        subplot(7,1,1)
        hold on
        plot(t,Vm,'b','LineWidth', 1.2)
        hold on
        plot(t(TF),Vm(TF),'*r')
        plot(t(TB),Vm(TB),'*g')

        %   plot(t(seuil_th), Vm(seuil_th),'sb','MarkerSize',10)


        hold off
        % plot(t(index_mdp(end)),Vm(index_mdp(end)),'r*','LineWidth', 1.2)
        ylabel('Vm (mV)')
        subplot(7,1, 2)
        hold on
        plot(t,ICaL,'b','LineWidth', 1.2)
        ylabel('ICaL (pA/pF)')
        subplot(7,1,3)

        plot(t,ICaT,'b','LineWidth', 1.2)
        ylabel('ICaT (pA/pF)')
        subplot(7,1, 4)
        hold on
        plot(t,ICaD,'b','LineWidth', 1.2)
        ylabel('ICaD (pA/pF)')

        subplot(7,1, 5)
        hold on
        plot(t,If,'b','LineWidth', 1.2)
        ylabel('If (pA/pF)')

        subplot(7,1, 6)
        hold on
        plot(t, I_Gap_AV,'b','LineWidth', 1.2)
        ylabel('I_Gap_AV (pA/pF)')


        subplot(7,1, 7)
        hold on
        plot(t, Istim,'b','LineWidth', 1.2)
        ylabel('IStim (nA/pF)')
        sgtitle([' PA Myo Gap = ',num2str(G_gap) ,' StimMac ',num2str(stim_flagM),...
            ' StimMyo ',num2str( stim_flag ),' Stim Frequency ',num2str(i_stim_frequency )])
        %
        % hold on
        saveas(gcf,[dir_result,'PA_Ggap=' num2str(G_gap) '_Ko=' num2str(K_e_Phi) '.jpg'],'jpg');
        saveas(gcf,[dir_result,'PA_Ggap=' num2str(G_gap) '_Ko=' num2str(K_e_Phi) '.fig'],'fig');
        close

        figure('units','normalized','outerposition',[0 0 1 1])

        subplot(8,1,1)
        hold on
        plot(t,VmMphi,'b','LineWidth', 1.2)
        hold off
        % plot(t(index_mdp(end)),Vm(index_mdp(end)),'r*','LineWidth', 1.2)
        ylabel('Vm Mphi (mV)')

        subplot(8,1, 2)
        hold on
        plot(t, I_Gap_Phi,'b','LineWidth', 1.2)
        ylabel('Igap (pA/pF)')

        subplot(8,1,3)
        plot(t,  Ib_Phi,'b','LineWidth', 1.2)
        ylabel('I K Back (pA/pF)')

        subplot(8,1, 4)
        hold on
        plot(t, IK1_Phi,'b','LineWidth', 1.2)
        ylabel('IK1_MPhi (pA/pF)')

        subplot(8,1, 5)
        hold on
        plot(t,IKur_Phi,'b','LineWidth', 1.2)
        ylabel('IKur_MPhi (pA/pF)')





        subplot(8,1, 6)
        hold on
        plot(t,Ishk_Phi,'b','LineWidth', 1.2)
        ylabel('Ishk_MPhi (pA/pF)')


        subplot(8,1, 7)
        hold on
        plot(t,(Ishk_Phi+IK1_Phi+Ib_Phi+IKur_Phi +I_Gap_Phi),'b','LineWidth', 1.2)
        ylabel('current net (nA/pF)')

        subplot(8,1, 8)
        hold on
        plot(t,i_diff_Phi,'b','LineWidth', 1.2)
        ylabel('Stim Mphi (nA/pF)')
        sgtitle([' PA Macro Gap = ',num2str(G_gap) ,' StimMac ',num2str(stim_flagM),' StimMyo ',...
            num2str( stim_flag ),' Stim Frequency ',num2str(i_stim_frequency )])
        saveas(gcf,[dir_result,'PAMacro_Ggap=' num2str(G_gap) '_Ko=' num2str(K_e_Phi) '.jpg'],'jpg');
        saveas(gcf,[dir_result,'PAMacro_Ggap=' num2str(G_gap) '_Ko=' num2str(K_e_Phi) '.fig'],'fig');
        close






        G_gap=G_gap+gapstep;


    end


    if resultat_sauv_affiche==1
        %[G_gap Choix_Type_Phi K_e_Phi  (sum(TF)/60) mean(diff(t(TF))) mean(Vm(TF)) mean(Vm(TB))];
        % ResultSimT=table(ResultVm(:,1),ResultVm(:,2),ResultVm(:,3),ResultVm(:,4),ResultVm(:,5),ResultVm(:,6),ResultVm(:,7),ResultVm(:,8),ResultVm(:,9),ResultVm(:,10));
        % ResultSimT.Properties.VariableNames=["G_Gap nS/pF","G_Gap nS","Type Macro","Concentration K ext","Rate (PA/min)","CL","PA maxMoy(mV)","MDP moy (mV)" "I gap moy (pA/pF)" "Amplitude moy PA"];
        % NomRes=[dir_result,'MacrophageMyo_Resultats','_type',num2str(Choix_Type_Phi),'.csv'];
        % writetable( ResultSimT,NomRes,'Delimiter','tab')


        %%
        figure('units','normalized','outerposition',[0 0 1 1])
        hold on
        plot3(ResultVm(:,1),ResultVm(:,4),ResultVm(:,7),liste_Sym{i_ke},...
            'Color',linecolors(i_ke,:),'MarkerSize',10,'MarkerFaceColor',linecolors(i_ke,:))
        set(gca,'FontSize',5)
        set(gca,'LineWidth',5)
        xlabel('Conductance Gap nS/pF','FontSize',20,'FontWeight','bold')
        ylabel('PA/min','FontSize',20,'FontWeight','bold')
        zlabel('Potentiel min Moyen (mV)' ,'FontSize',20,'FontWeight','bold')
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',20,'FontWeight','bold')
        saveas(gcf,[dir_result,'PA_Ggap_Vm=' num2str(G_gap) '_Ko=' num2str(K_e_Phi) '_' num2str(Choix_Type_Phi) '.jpg'],'jpg');

        view([-41 55])
        saveas(gcf,[dir_result,'PA_Ggap_Vm3D=' num2str(G_gap) '_Ko=' num2str(K_e_Phi) '_' num2str(Choix_Type_Phi) '.jpg'],'jpg');
        close
        figure('units','normalized','outerposition',[0 0 1 1])
        hold on
        plot(ResultVm(:,1),ResultVm(:,7),liste_Sym{i_ke},...
            'Color',linecolors(i_ke,:),'MarkerSize',10,'MarkerFaceColor',linecolors(i_ke,:))
        set(gca,'FontSize',5)
        set(gca,'LineWidth',5)
        xlabel('Conductance Gap nS/pF','FontSize',20,'FontWeight','bold')

        ylabel('Potentiel min Moyen (mV)' ,'FontSize',20,'FontWeight','bold')
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',20,'FontWeight','bold')
        saveas(gcf,[dir_result,'PA_Ggap_MDP=' num2str(G_gap) '_Ko=' num2str(K_e_Phi) '_' num2str(Choix_Type_Phi) '.jpg'],'jpg');
        %legend({'Ko = 5.2 mM','Ko=10 mM','Ko=50 mM','Ko=120 mM'})
        close
    end
end
