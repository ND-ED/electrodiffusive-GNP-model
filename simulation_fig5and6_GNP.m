F=96485;s_r_v=2/(5*10^(-6));gamma=(s_r_v/F)*10^(-2);
c_K_t=100;c_Cl_t=145;c_Na_t=155;tao=1000;belta=1;
P_KL=s_r_v*2*10^(-8);P_NaL=s_r_v*1*10^(-9);P_ClL=s_r_v*4*10^(-9);

pc=((30/(4*pi*(2.5)^2))*26.64/(9648500))*(10^7);
c_Na_i=XX3(5);c_Na_o=155-c_Na_i;c_K_o=XX3(4);c_K_i=100-c_K_o;
r=(c_K_i-c_K_o)/(c_Na_o-c_Na_i);c_Gl_o=c_K_o+r*c_Na_o;c_Gl_i=c_K_i+r*c_Na_i;
Glbart=c_Gl_o;P_Kgl=pc/Glbart;P_Nagl=r*P_Kgl;

P_gluK=s_r_v*P_Kgl*pi*(10^(-8));P_gluNa=s_r_v*P_Nagl*pi*(10^(-8));glu=0;

I_pump=0.048;g_KV=20;g_NaV=30;
l_t=1500000;dt=0.1;%fig. 5g: l_t=50000; %fig. 5jk and fig. 6: l_t=1500000;
stepc=1;t_tsee=stepc:stepc:l_t;l_tc=length(t_tsee);
X_inp=XX3;dA=-110;

X_reco1=zeros(l_tc,7);X_RTM=zeros(1,6);A_reco=zeros(l_tc,1);Y_RTM=zeros(1,11);
G_RTM(2)=P_gluK;G_RTM(3)=P_gluNa;

for i2=1:1

    X_RTM(:)=X_inp;
    
    for i1=1:l_t

        tt1=0.2*l_t;tt2=l_t;
        if (i1>tt1)&&(i1<tt2)
            glu=0.515;%fig. 5g 0.195 fig. 5j and fig. 6 0.515 fig. 5k 0.525
        else
            glu=0;
        end
        G_RTM(1)=glu;

        Y_RTM(:)=[P_KL,P_NaL,P_ClL,c_K_t,c_Na_t,c_Cl_t,tao,dA,gamma,g_KV,g_NaV];
        k1=RTM_GNP_b(X_RTM,Y_RTM,Z_RTM,G_RTM);k2=RTM_GNP_b(X_RTM+(dt/2)*k1,Y_RTM,Z_RTM,G_RTM);
        k3=RTM_GNP_b(X_RTM+(dt/2)*k2,Y_RTM,Z_RTM,G_RTM);k4=RTM_GNP_b(X_RTM+dt*k3,Y_RTM,Z_RTM,G_RTM);

        X_RTM=X_RTM+(dt/6)*(k1+2*k2+2*k3+k4);

        for i4=1:3
            if X_RTM(i4)>1
                X_RTM(i4)=1;
            elseif X_RTM(i4)<0
                X_RTM(i4)=0;
            end
        end

        if rem(i1,stepc)==0
            ii1=i1/stepc;
            X_reco1(ii1,2:7)=X_RTM;
            A_reco(ii1,:)=dA;
        end

    end
    X_reco1(:,1)=(1000/gamma)*(X_reco1(:,6)-X_reco1(:,5)-X_reco1(:,7)+(c_K_t-c_Na_t+c_Cl_t+dA)/2);

end