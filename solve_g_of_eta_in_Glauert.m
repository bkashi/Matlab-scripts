% In Glauert 1956 f'=2gg', where g is given implicitly as g(eta).
% Here we plot f'(eta) from Glauert wall jet (1956).
% First we find g(eta) numerically,
% then differentiate to obtain g'(eta),
% then plot f'(eta) for several grids sizes,
% Finally, export eta and f'(eta) to a csv file.

clear variables
close all

% --------------------------------- %
% grid ind. test for numerical sol. %
% --------------------------------- %
fprintf('\nN \t \t u \t\t\t eta\n')
fprintf('------ \t ----- \t\t -----\n')

figure('position',[200,200,300,250])
hold on; box on

symbols=['b- '; 'g- ' ;'--r'];
widths=[2 1.5 1.5];
N=[100 1000 10000];
for i=1:length(N) % discretization numer of points
  n=N(i);
  g=linspace(0,1,n);
  eta=log(sqrt(g.^2+g+1)./(1-g))+sqrt(3)*atan(sqrt(3)*g./(2+g));
  g_=[diff(g)./diff(eta) 0];
  
  f_=2*g.*g_;
  plot(eta,f_,symbols(i,:),'linewidth',widths(i))

  [a,b]=max(f_);
  fprintf('%d \t %1.4f \t %1.4f \n', n, a, eta(b))
end
xlabel('\eta','fontsize',10)
ylabel("f '",'fontsize',10)
legend('n=100','     1,000','     10,000')

disp('')

% write files
%dlmwrite('numerical_Glauert_fd.txt',[eta(1:10:end)' f_(1:10:end)'],' ')
%f__=diff(f_)./diff(eta);
%dlmwrite('numerical_Glauert_fdd.txt',[eta(1:10:end-1)' f__(1:10:end)'],' ')
















if (1==0)
% ---------------------- %
% compare approximations %
% ---------------------- %
figure; hold on; box on
xlim([0 10])
ylim([0 0.35])
xlabel('\eta','fontsize',14)
ylabel("f '",'fontsize',15)
ax = gca;
ax.XAxisLocation = 'origin';

% numerical
plot(eta,f_,'-g','linewidth',1.4)

% Mo et al.
plot(eta,2*0.125.*eta.*exp(-0.125*eta.^2),'b','linewidth',1.8)

% Schwarz and Caswell
eta_m=0.4831*eta;
plot(eta,0.315*(0.4831*eta).*(1-0.2334*eta.^2),'m')

% Bouremel exact (sampled from fig. 3)
Bouremel_x=[0.05627	0.003139	0.117	0.1976	0.1705	0.1519	0.2933	0.2635	0.3937	0.3472	0.496	0.4402	0.5983	0.5518	0.6987	0.6448	0.7935	0.7377	0.8989	0.8493	0.9981	0.9516	1.097	1.054	1.195	1.156	1.298	1.249	1.4	1.5	1.597	1.733	1.881	2.086	2.263	2.384	2.504	2.579	2.673	2.743	2.802	2.89	2.989	3.088	3.146	3.197	3.286	3.386	3.491	3.591	3.658	3.723	3.799	3.899	3.994	4.095	4.19	4.29	4.401	4.504	4.606	4.736	4.857	4.978	5.136	5.341	5.527	5.713	5.917	6.122	6.326	6.531	6.735	6.94	7.176	7.572	7.971	8.167	8.372	8.577	8.781	8.986	9.19	9.395	9.599	9.804	9.916];
Bouremel_y=[0.01563	0.004989	0.02535	0.04821	0.04136	0.03418	0.06615	0.05952	0.08791	0.07812	0.1102	0.09865	0.1315	0.1227	0.1523	0.1424	0.173	0.1618	0.1935	0.1844	0.2116	0.2042	0.2298	0.2218	0.2463	0.2398	0.2613	0.2545	0.2747	0.2863	0.2964	0.307	0.3121	0.3137	0.3097	0.3022	0.2943	0.287	0.2772	0.2718	0.2647	0.2545	0.2425	0.2308	0.2226	0.2172	0.2057	0.1933	0.1809	0.1692	0.1617	0.1544	0.1452	0.1353	0.1245	0.1158	0.1069	0.09861	0.08982	0.08131	0.0754	0.06624	0.06098	0.05355	0.04667	0.03811	0.03255	0.02683	0.02205	0.01803	0.01469	0.01204	0.00963	0.007544	0.005957	0.00401	0.002459	0.001976	0.001787	0.001365	0.0009984	0.0005034	0	0	0	0	0];
plot(Bouremel_x,Bouremel_y,'--k','linewidth',1.2) % sampled from fig. 3 in Bouremel

% Haustein
 plot(eta,2*erf(eta/pi).*(2*exp(-eta.^2/pi^2)/pi^(3/2)),'r:','linewidth',1.5)
 
legend('numerical','Mo et al.','Schwarz and Caswell','Bouremel (sampled)', 'new')

% ------------------- %
% compare derivatives %
% ------------------- %
figure; hold on; box on
xlim([0 10])
ylim([-0.15 0.3])
xlabel('\eta','fontsize',14)
ylabel("f '",'fontsize',15)
ax = gca;
ax.XAxisLocation = 'origin';

% numerical
plot(eta(2:end),diff(f_)./diff(eta),'-g','linewidth',1.4)

% Mo et al.
plot(eta,0.25*exp(-0.125*eta.^2)-0.0625.*eta.^2.*exp(-0.125.*eta.^2),'b','linewidth',1.8)

% Schwarz and Caswell
plot(eta,-.10655*eta.^2+.15218,'m')

% Haustein
plot(eta,8*(exp(-eta.^2/pi^2)).^2/pi^3-8*erf(eta/pi).*eta.*exp(-eta.^2/pi^2)/pi^(7/2),'r:','linewidth',1.5)
 
legend('numerical','Mo et al.','Schwarz and Caswell','new')
end
