clear all; close all

tau=1;%0.001*sqrt(1000/101325); % =R_inf/sqrt(rho/P_inf)

[t,y] = ode45(@eqs2,[0 0.0025],[0; 0.00101; 104369.5]);

plot(t/tau,y(:,2))

xlabel('\it t \rm / (\it R_{\infty}\rm /\surd\rho/\it g_0P_{\infty})\rm [-]')
ylabel('\it R / R_{\infty} \rm [-]')

%xlim([0 50])
%ylim([0.9 1.1])