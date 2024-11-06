F=96485;s_r_v=2/(5*10^(-6));
gamma=(s_r_v/F)*10^(-2);
g_Na=30;g_K=20;
c_K_t=100;c_Cl_t=145;c_Na_t=155;
tao=1000;belta=1;dA=-110;


X=XX3;
m=X(1);h=X(2);n=X(3);c_K_o=X(4);c_Na_i=X(5);c_Cl_i=X(6);
c_K_i=c_K_t-c_K_o;c_Na_o=c_Na_t-c_Na_i;c_Cl_o=c_Cl_t-c_Cl_i;
V=(1000/gamma)*(c_Na_i-c_K_o-c_Cl_i+(c_K_t-c_Na_t+c_Cl_t+dA)/2);
E_K=26.64*log(c_K_o/c_K_i);E_Na=26.64*log(c_Na_o/c_Na_i);E_Cl=26.64*log(c_Cl_i/c_Cl_o);
X_inp=[V,m,h,n];

Kbar0=(c_K_i-c_K_o*exp(-V/26.64))*(1/(1-exp(-V/26.64)))*(V/(V-26.64*log(c_K_o/c_K_i)));
Nabar0=(c_Na_i-c_Na_o*exp(-V/26.64))*(1/(1-exp(-V/26.64)))*(V/(V-26.64*log(c_Na_o/c_Na_i)));
Clbar0=(c_Cl_i-c_Cl_o*exp(V/26.64))*(1/(1-exp(V/26.64)))*(V/(V-26.64*log(c_Cl_i/c_Cl_o)));
P_KL=2*10^(-8);P_NaL=1*10^(-9);P_ClL=4*10^(-9);
I_pump=0.048;G_glu=0.12;E_glu=0;

l_t=200000;%fig. 4: l_t=1000;fig. 5 a,d: l_t=200000;
dt=0.1;
stepc=1;t_tsee=stepc:stepc:l_t;l_tc=length(t_tsee);

X_recol=zeros(l_tc,4);X_RTM=zeros(1,4);Y_RTM=zeros(1,9);L_RTM=zeros(1,9);

for i2=1:1

    X_RTM(:)=X_inp;

    if i2==1
        feq_gaba=80;
    else
        feq_gaba=25;
    end
    
    for i1=1:l_t
        
        tt1=0.2*l_t;tt2=l_t;
        if i1>tt2
            Glu=0;
        elseif i1>tt1
            Glu=0.195; 
        else
            Glu=0;
        end


        Y_RTM(:)=[g_K,g_Na,E_K,E_Na,gamma,I_pump,Glu*G_glu,E_glu,F];
        L_RTM(:)=[c_K_i,c_K_o,c_Na_i,c_Na_o,c_Cl_i,c_Cl_o,P_KL,P_NaL,P_ClL];

        k1=RTM_fGNP(X_RTM,Y_RTM,L_RTM);k2=RTM_fGNP(X_RTM+(dt/2)*k1,Y_RTM,L_RTM);
        k3=RTM_fGNP(X_RTM+(dt/2)*k2,Y_RTM,L_RTM);k4=RTM_fGNP(X_RTM+dt*k3,Y_RTM,L_RTM);

        for i4=2:4
            if X_RTM(i4)>1
                X_RTM(i4)=1;
            elseif X_RTM(i4)<0
                X_RTM(i4)=0;
            end
        end
        X_RTM=X_RTM+(dt/6)*(k1+2*k2+2*k3+k4);

        if rem(i1,stepc)==0
            ii1=i1/stepc;
            X_recol(ii1,:)=X_RTM;
        end       
    end

end
