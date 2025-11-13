clc
clear
close all
freq_series=[0,350,450,550,720];
typemodel=[0,2];
for freq=1:1%2:length(freq_series)
    for tm=1:length(typemodel)
        freqbpm=freq_series(freq);
        Cm_Phi=18.32;
        Choix_Type_Phi=typemodel(tm);
        StimMac=1;
        stmac=0;
        stmyo=0;
        dGap=0.05;
        K_e_Phi=5.4;
        nsec = 20; %s
        CL1= nsec;
        perc_Amp=0.65;
        fenetre_detect=0.02;
        Capa=18.32;
        stimtype='_Spont_';
        %stimtype='_StimMac_';
        %stimtype='_StimMyo_';
        methodedetect_seuil=0;
        minPerioddetect=0.06;
        seuildvdt=0.19;clc
        
        session=150;
        %seuildvdt=2;
        Nom_fichier_Entree=['F:\ResultatsSimulationMacrophages\Resultats03112025_typeMPhi',...
            num2str(Choix_Type_Phi),...
            'GgapStep',num2str(dGap),'StimMacro',num2str(stmac),'StimMyo',...
            num2str(stmyo),'perStim',num2str(freqbpm),'AmpSt600session',num2str(session),...
            'capaM',num2str(Capa),'\'];



        dt=0.00025;
        PA=[];time_n = 0:dt:nsec; time_n=time_n';

        [fname,chemin]=uigetfile([Nom_fichier_Entree,'*.csv'],'MultiSelect','on');%loads file



        if isequal([fname chemin],[0,0])
            return
        else


            session=size(fname,2);
            Ggap=(0:session).*dGap;
            PAMatrix=zeros(length(time_n),size(fname,2));
            dVmMatrix=zeros(length(time_n),size(fname,2));
            for i=1:size(fname,2)

                datname=[chemin,fname{i}];
                PACur=readtable(datname);%,'VariableNamingRule','preserve');
                time=PACur.time;
                PA=PACur.Vm;

                PA_int = spline(time,PA,time_n);

                PAMatrix(:,i)=PA_int;

                %dVmMatrix(:,i)  = [0;diff(PA_int)./diff(time_n)];
                dVmMatrix(:,i)  = [0;diff(PA_int)];

            end

            path=pwd();
            %dir_result=[ path,'\ResultatsAnalyse_typeMPhi',num2str( Choix_Type_Phi),'StimMacro',num2str( stim_flagM),...
            %   'StimMyo',num2str( stim_flag),'gapstep005nSsession',num2str(session),'capaM',num2str(Cm_Phi),'\'];
            dir_result=[ 'F:\ResultatsSimulationMacrophages\ResultatsAnalyse03112025_typeMPhiv',num2str( Choix_Type_Phi),...
                stimtype,num2str(StimMac),'_freq',num2str( freqbpm),'_',strtok(fname{session},'.'),'\'];
            test_dir=exist(dir_result,'dir');
            if test_dir==0
                mkdir(dir_result)
            end


            figure('units','normalized','outerposition',[0 0 1 1])


            imagesc(time_n,Ggap,PAMatrix')

            saveas(gcf,[dir_result,'PAMATRIX_', stimtype,num2str(StimMac),...
                '_type', num2str(Choix_Type_Phi),'BPM',num2str(freqbpm) ,'.jpg'],'jpg');
            close
            %%

            Result_Vm=zeros(1,10);
            %********************************************************************************************


            %*********************************************************************************************

            %PAR_PA=[];

            PAR_PA= table('Size',[10000 22],'VariableTypes',{'double','double','double','double','double','double', ...
                'double','double','double','double','double','double','double',...
                'double','double','double','double','double','double','double','double','double'},...
                'VariableNames',...
                {'NoEpisode','NoPA','ASeuil','TSeuil','AMax','TMax','MDP', 'TMDP'...
                'A10Monte','T10Monte','A50Monte','T50Monte','T90Monte','A90Monte',...
                'A10Repol', 'T10Repol','A50Repol','T50Repol','A90Repol','T90Repol',...
                'DureePA50','DureePA90'}) ;




            count=1;


            %%


            for sim=1:session
                Vm=PAMatrix(:,sim);
                Vmdetect=Vm+100;
                maxAmpdetect=max(Vmdetect)-max(Vmdetect)* perc_Amp;

                dVm  =dVmMatrix(:,sim);% [0;diff(Vm)./diff(time_n)];
                maxAmpdvdetect=max( dVm)-max( dVm)*perc_Amp;
                TF = islocalmax(Vmdetect,'MinProminence', maxAmpdetect,'MinSeparation', minPerioddetect,'SamplePoints',time_n);
                TB= islocalmin(Vmdetect,'MinProminence', maxAmpdetect,'MinSeparation', minPerioddetect,'SamplePoints',time_n);

                %TFd=islocalmax(dVm,'MinProminence',maxAmpdvdetect,'MinSeparation', minPerioddetect,'SamplePoints',time_n);
                indexmaxPA=find(TF);
                indexMDPPA=find( TB);
                % indexMaxdvdt=find(TFd);
                CL=mean(diff(time_n(TF)));
                indexmaxPA=indexmaxPA(2:end-1);

                %indexMaxdvdt=indexMaxdvdt(2:end-1);

                %indexMaxdvdt=indexMaxdvdt(2:end-1);

                for imax=1:length(indexmaxPA)
                    index=find(indexMDPPA> indexmaxPA(imax),1,'first');
                    if ~isempty(index)
                        indexMDPPA1(imax)=  indexMDPPA(index);


                    end
                end
                indexMDPPA= indexMDPPA1;


                %        indexmax=find(TFd);
                if  (~isempty(TF)) && sum(TF)>2
                    for i=1:length(indexmaxPA)
                        tpeak=indexmaxPA(i);
                        tmdp=indexMDPPA(i);


                        tss=tpeak-round(fenetre_detect/dt);


                        Amp_PAD=Vm(tpeak)-Vm(tmdp);
                        SAD90=Vm(tpeak)-Amp_PAD*0.9;
                        if  methodedetect_seuil==0

                            %[~,peakdvdt_t]=max(dVm(tss:tpeak));


                            for j=tss:-1:1

                                if  dVm(j)<seuildvdt
                                    tseuil=j;
                                    break;

                                end
                            end

                        elseif methodedetect_seuil==1

                            [~,peakdvdt_t]=max(dVm(tss:tpeak));



                            for j=(tpeak-20):-1:1

                                if (dVm(j)<dVm(j-1))
                                    tseuil=j;
                                    break;

                                end
                            end




                        elseif methodedetect_seuil==2

                            idown=find(Vm(tss:tpeak)>SAD90,1,"first");
                            tseuil=tpeak-idown;

                        end




                        Amp_PAM=Vm(tseuil)-Vm(tpeak);
                        SAM10=Vm(tseuil)+Amp_PAM*0.1;
                        SAM50=Vm(tseuil)+Amp_PAM*0.5;
                        SAM90=Vm(tseuil)+Amp_PAM*0.9;
                        SAD10=Vm(tpeak)-Amp_PAD*0.1;
                        SAD50=Vm(tpeak)-Amp_PAD*0.5;






                        PAR_PA.NoEpisode(count)=sim;
                        PAR_PA.NoPA(count)=i;
                        PAR_PA.ASeuil(count)=Vm(tseuil);
                        PAR_PA.TSeuil(count)=time_n(tseuil);
                        PAR_PA.AMax(count)=Vm(  tpeak);
                        PAR_PA.TMax(count)=time_n(  tpeak);
                        PAR_PA.MDP(count)=Vm(tmdp);

                        PAR_PA.TMDP(count)=time_n(tmdp);
                        AmpPA(count)=Vm(tpeak)-Vm(tmdp);




                        iup=find(Vm(tseuil:tpeak)>SAM10,1,"first");
                        idown=find(Vm( tpeak: tmdp)>SAD10,1,"last");
                        tup=tseuil+iup;
                        tdown= tpeak+idown;
                        PAR_PA.A10Monte(count)=Vm(tup);
                        PAR_PA.T10Monte(count)=time_n(tup);
                        PAR_PA.A10Repol(count)=Vm(tdown);
                        PAR_PA.T10Repol(count)=time_n(tdown);

                        iup=find(Vm(tseuil:  tpeak)>SAM50,1,"first");
                        idown=find(Vm( tpeak: tmdp)>SAD50,1,"last");
                        tup=tseuil+iup;
                        tdown= tpeak+idown;
                        PAR_PA.A50Monte(count)=Vm(tup);
                        PAR_PA.T50Monte(count)=time_n(tup);
                        PAR_PA.A50Repol(count)=Vm(tdown);
                        PAR_PA.T50Repol(count)=time_n(tdown);

                        iup=find(Vm(tseuil:tpeak)>SAM90,1,"first");
                        idown=find(Vm(tpeak: tmdp)>SAD90,1,"last");
                        tup=tseuil+iup;
                        tdown= tpeak+idown;
                        PAR_PA.A90Monte(count)=Vm(tup);
                        PAR_PA.T90Monte(count)=time_n(tup);
                        PAR_PA.A90Repol(count)=Vm(tdown);
                        PAR_PA.T90Repol(count)=time_n(tdown);
                        PAR_PA.DureePA50(count)= PAR_PA.T50Repol(count)-PAR_PA.T50Monte(count);
                        PAR_PA.DureePA90(count)=  PAR_PA.T90Repol(count)-PAR_PA.TSeuil(count);


                        count=count+1;
                    end



                end
                %   ResultSimT.Properties.VariableNames=["G_Gap nS/pF","G_Gap nS","Type Macro","Concent_ion K ext",...
                %"Rate (PA/min)","CL","PA maxMoy(mV)","MDP moy (mV)" "Amplitude moy PA" "APD90"];

                if  (~isempty(TF)) && sum(TF)>2%(~isempy(TB))
                    Result_Vm(sim,:)=[(Ggap(sim)/Cm_Phi) Ggap(sim) Choix_Type_Phi K_e_Phi  ...
                        (sum(TF)*60/CL1) mean(diff(time_n(TF))) mean(Vm(TF)) mean(Vm(TB))   ...
                        mean(AmpPA)  mean(PAR_PA.DureePA90( PAR_PA.NoEpisode==sim)) ];
                elseif sum(TF)<=2

                    Result_Vm(sim,:)=[(Ggap(sim)/Cm_Phi) Ggap(sim) Choix_Type_Phi K_e_Phi  0 0 mean(Vm(2000:end)) mean(Vm(2000:end))   0  0];
                    % else
                    %     ResultVm(sim,:)=[(G_gap/Cm_Phi) Choix_Type_Phi K_e_Phi  0 0 mean(Vm(2000:end)) 0 mean(I_Gap_Phi) 0];
                end

            end
            PAR_PA(count:end,:)=[];

            NomRes=[dir_result,'ParametresPA_Resultats','_type',num2str(Choix_Type_Phi),stimtype,'BPM',num2str( freqbpm),'.csv'];
            writetable( PAR_PA,NomRes,'Delimiter','tab')



        end


        %%
        for i=1:session
            num=i;
            tseuil=round(PAR_PA.TSeuil(PAR_PA.NoEpisode==num)./dt);
            tmax=round(  PAR_PA.TMax(PAR_PA.NoEpisode==num)./dt);
            t90=round(PAR_PA.T90Repol(PAR_PA.NoEpisode==num)./dt);
            tmdp=round(PAR_PA.TMDP(PAR_PA.NoEpisode==num)./dt);


            figure(i)
            sgtitle(['Ggap ',num2str(Ggap(i)),' ',stimtype,' ',' BPM ',   num2str(freqbpm)])
            plot(time_n,PAMatrix(:,num),'b')
            hold on
            plot(time_n(tseuil),PAMatrix(tseuil,num),'dk','MarkerFaceColor','k')
            plot(time_n(t90),PAMatrix(t90,num),'>r')
            plot(time_n(tmax),PAMatrix(tmax,num),'*r')
            plot(time_n(tmdp),PAMatrix(tmdp,num),'or')

            yyaxis("right")
            plot(time_n,dVmMatrix(:,num),'r')
            hold on

            %pause
            saveas(gcf,[dir_result,'PA_','G_gap',num2str(Ggap(i)),'_type', num2str(Choix_Type_Phi) ,'.fig'],'fig');
            close
        end




        ResultSimT=table(Result_Vm(:,1),Result_Vm(:,2),Result_Vm(:,3),Result_Vm(:,4),Result_Vm(:,5),...
            Result_Vm(:,6),Result_Vm(:,7),Result_Vm(:,8),Result_Vm(:,9) ,Result_Vm(:,10));
        ResultSimT.Properties.VariableNames=["G_Gap nS/pF","G_Gap nS","Type Macro","Concent_ion K ext",...
            "Rate (PA/min)","CL","PA maxMoy(mV)","MDP moy (mV)" "Amplitude moy PA" "APD90"];

        NomRes=[dir_result,'MacrophageMyo_Resultats_type',num2str(Choix_Type_Phi),stimtype,'BPM',num2str( freqbpm),'.csv'];
        writetable( ResultSimT,NomRes,'Delimiter','tab')
    end
end