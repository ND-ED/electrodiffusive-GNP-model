function [XX] = RTM_fGNP(X,Y,L)
XX=zeros(1,4);

m=X(2);h=X(3);n=X(4);V=X(1);

g_K=Y(1);g_Na=Y(2);E_K=Y(3);E_Na=Y(4);gamma=Y(5);
I_pump=Y(6);G_glu=Y(7);E_glu=Y(8);F=Y(9);
c_K_i=L(1);c_K_o=L(2);c_Na_i=L(3);c_Na_o=L(4);c_Cl_i=L(5);c_Cl_o=L(6);
P_KL=L(7);P_NaL=L(8);P_ClL=L(9);

Kbart=(c_K_i-c_K_o*exp(-V/26.64))*(1/(1-exp(-V/26.64)))*(V/(V-26.64*log(c_K_o/c_K_i)));
Nabart=(c_Na_i-c_Na_o*exp(-V/26.64))*(1/(1-exp(-V/26.64)))*(V/(V-26.64*log(c_Na_o/c_Na_i)));
Clbart=(c_Cl_i-c_Cl_o*exp(V/26.64))*(1/(1-exp(V/26.64)))*(V/(V-26.64*log(c_Cl_i/c_Cl_o)));

G_KL=F*(100/26.64)*P_KL*Kbart;G_NaL=F*(100/26.64)*P_NaL*Nabart;G_ClL=F*(100/26.64)*P_ClL*Clbart;
G_L=G_KL+G_NaL+G_ClL;E_L=(G_KL*E_K+G_NaL*E_Na+G_ClL*(26.64*log(c_Cl_i/c_Cl_o)))/(G_KL+G_NaL+G_ClL);

XX(1)=-(g_Na*(m^3)*h)*(V-E_Na)-(g_K*(n^4))*(V-E_K)-G_glu*(V-E_glu)-G_L*(V-E_L)-I_pump/gamma;
XX(2)=((0.32*(V+54))/(1-exp(-(V+54)/4)))*(1-m)-((0.28*(V+27))/(exp((V+27)/5)-1))*m;
XX(3)=(0.128*exp(-(V+50)/18))*(1-h)-(4/(1+exp(-(V+27)/5)))*h;
XX(4)=((0.032*(V+52))/(1-exp(-(V+52)/5)))*(1-n)-(0.5*exp(-(V+57)/40))*n;

%if XX(2)/10>(1-m)
%    XX(2)=(0.9999-m)*10;
%end

%if XX(3)/10>(1-h)
%    XX(3)=(0.9999-h)*10;
%end

%if XX(4)/10>(1-n)
%    XX(4)=(0.9999-n)*10;
%end

end