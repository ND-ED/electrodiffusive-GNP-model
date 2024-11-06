function [XX] = RTM_chh(X,Y)
XX=zeros(1,4);

m=X(2);h=X(3);n=X(4);V=X(1);

g_K=Y(1);g_Na=Y(2);E_K=Y(3);E_Na=Y(4);gamma=Y(5);
I_pump=Y(6);G_glu=Y(7);E_glu=Y(8);G_L=Y(9);E_L=Y(10);

XX(1)=-(g_Na*(m^3)*h)*(V-E_Na)-(g_K*(n^4))*(V-E_K)-G_glu*(V-E_glu)-G_L*(V-E_L)-I_pump/gamma;
XX(2)=((0.32*(V+54))/(1-exp(-(V+54)/4)))*(1-m)-((0.28*(V+27))/(exp((V+27)/5)-1))*m;
XX(3)=(0.128*exp(-(V+50)/18))*(1-h)-(4/(1+exp(-(V+27)/5)))*h;
XX(4)=((0.032*(V+52))/(1-exp(-(V+52)/5)))*(1-n)-(0.5*exp(-(V+57)/40))*n;

end