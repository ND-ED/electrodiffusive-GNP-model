function [XX] = RTM_GNP_b(X,Y,G)
%UNTITLED2 此处提供此函数的摘要
%   此处提供详细说明
XX=zeros(1,3);

m=X(1);h=X(2);n=X(3);c_K_o=X(4);c_Na_i=X(5);c_Cl_i=X(6);

P_KL=Y(1);P_NaL=Y(2);P_ClL=Y(3);c_K_t=Y(4);c_Na_t=Y(5);c_Cl_t=Y(6);
tao=Y(7);dA=Y(8);gamma=Y(9);g_KV=Y(10);g_NaV=Y(11);

c_K_i=c_K_t-c_K_o;c_Na_o=c_Na_t-c_Na_i;c_Cl_o=c_Cl_t-c_Cl_i;

V=(1000/gamma)*(c_Na_i-c_K_o-c_Cl_i+(c_K_t-c_Na_t+c_Cl_t+dA)/2);


dnavt=-gamma*g_NaV*(m^3)*h*(V-26.64*log(c_Na_o/c_Na_i));dkvt=gamma*g_KV*(n^4)*(V-26.64*log(c_K_o/c_K_i));

P_gluK=G(1)*G(2);P_gluNa=G(1)*G(3);

dkt=(P_KL+P_gluK)*(V/26.64)*(c_K_i-c_K_o*exp(-V/26.64))/(1-exp(-V/26.64));

if V==0
    dkt=(P_KL+P_gluK)*(c_K_i-c_K_o);
end

dnat=-(P_NaL+P_gluNa)*(V/26.64)*(c_Na_i-c_Na_o*exp(-V/26.64))/(1-exp(-V/26.64));
if V==0
    dnat=-(P_NaL+P_gluNa)*(c_Na_i-c_Na_o);
end

dclt=P_ClL*(V/26.64)*(c_Cl_i-c_Cl_o*exp(V/26.64))/(1-exp(V/26.64));
if V==0
    dclt=-P_ClL*(c_Cl_i-c_Cl_o);
end

I_pump=(0.299510717130278/(1+exp((25-c_Na_i)/3)))*(1/(1+exp(3.5-c_K_o)));

XX(1)=((0.32*(V+54))/(1-exp(-(V+54)/4)))*(1-m)-((0.28*(V+27))/(exp((V+27)/5)-1))*m;
XX(2)=(0.128*exp(-(V+50)/18))*(1-h)-(4/(1+exp(-(V+27)/5)))*h;
XX(3)=((0.032*(V+52))/(1-exp(-(V+52)/5)))*(1-n)-(0.5*exp(-(V+57)/40))*n;
XX(4)=(dkt+dkvt-2*I_pump)/tao;
XX(5)=(dnat+dnavt-3*I_pump)/tao;
XX(6)=(dclt)/tao;

end