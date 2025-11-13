function [start_cal , End_cal]=Recher_Start_Stop(data,pidx,seuil_c,Fenetre_detection)


                    

                    start_cal=0;End_cal=0;

               
                    



                        s_f=(pidx-Fenetre_detection);
                        if s_f<1
                            s_f=1;
                        end
                        e_f= (pidx);%
                        if e_f>size(data,1)
                            e_f=size(data,1);
                        end

                       % [~,I_max_dvdt]=max(data(s_f:e_f));
                        [~,I_min_dvdt]=min(data(s_f:e_f));
                       
                       % ref_rech_max=pidx-I_max_dvdt;
                        ref_rech_min=pidx+I_min_dvdt;
                   
                        for k=s_f:-1:1%

                            if data(k)<seuil_c
                                start_cal=k;
                                break;

                            end
                        end


                        for k=ref_rech_min:1:size(data,1)%

                            if data(k)<seuil_c
                                End_cal=k;
                                break;

                            end

                        end
                    end
