clear variables
close all

Re=375; % Reynolds in experiment
figure('position',[150 40 600 400]); hold on
zone=3;
switch zone
    case 0 % stagnation - uniform and developed
        xlim([-0.5 1.2]);
        ylim([0 0.15]);
        Glauert_R=        [];
        Schlichting_R=    [0.05 0.1 0.15 0.25 0.4];
        Schlichting_col=  ['r' 'g' 'b' 'm' 'k' 'c'];
        sim_dev_R=        [0.05 0.1 0.15 0.25 0.4];
        sim_dev_col=      ['r' 'g' 'b' 'm' 'k' 'c'];
        experimental_R=   [];
        hand1=plot([0 0],[-1 -2],'k-'); % legend simulation
        hand2=plot([0 0],[-1 -2],'k--'); % legend Schlichting
        hands=[hand1(1) hand2(1)];
        plot([0 0],[0 0.2],'k.-.','linewidth',1.2) % centerline
        S=3; % data reduction
        H=0.3; % position of annotations
        V=0.05;
        
    case 1 % stagnation
        xlim([0 1.5]);
        ylim([0 0.15]);
        Glauert_R=        [];
        Schlichting_R=    [0.25 0.4 0.7];
        Schlichting_col=  ['r' 'g' 'b' 'm' 'k' 'c'];
        sim_dev_R=        [0.25 0.4 0.7];
        sim_dev_col=      ['r' 'g' 'b' 'm' 'k' 'c'];
        experimental_R=   [0.25 0.4 0.7];
        experimental_col= ['r' 'g' 'b'];
        hand1=plot([0 0],[-1 -2],'k-'); % legend numerical
        hand2=plot([0 0],[-1 -2],'ko'); % legend experimental
        hand3=plot([0 0],[-1 -2],'k--'); % legend Schlichting
        hands=[hand1(1) hand2(1) hand3(1)];
        S=1; % data reduction
        H=0.3;
        V=0.05;
        
    case 2 % transition
        xlim([0 1.5]);
        ylim([0 0.5]);
        Glauert_R=        [1.5];
        Glauert_col=      ['m'];
        sim_dev_R=        [0.7 0.9 1.15 1.5];
        sim_dev_col=      ['r' 'g' 'b' 'm'];
        experimental_R=   [0.7 0.9 1.15];
        experimental_col= ['r' 'g' 'b'];
        Schlichting_R=    [];
        Schlichting_col=  [];
        hand1=plot([0 0],[-1 -2],'k-'); % legend numerical
        hand2=plot([0 0],[-1 -2],'ko'); % legend experimental
        hand3=plot([0 0],[-1 -2],'k--'); % legend Glauert
        hands=[hand1(1) hand2(1) hand3(1)];
        S=1;
        H=0.3;
        V=0.3;
        
    case 3 % far wall
        Glauert_R=        [1  2  3.6  5];
        Glauert_col=      ['r' 'g' 'b' 'm' 'k'];
        Schlichting_R=    []; 
        sim_dev_R=        [1  2  3.6  5];
        sim_dev_col=      ['r' 'g' 'b' 'm' 'k'];
        experimental_R=   [1  2  3.6];
        experimental_col= ['r' 'g' 'b'];
        xlim([0 1.2]);
        ylim([0 0.4]);
        hand1=plot([0 0],[-1 -2],'k-'); % legend numerical
        hand2=plot([0 0],[-1 -2],'ko'); % legend experimental
        hand3=plot([0 0],[-1 -2],'k--'); % legend Glauert/Schlichting
        hands=[hand1(1) hand2(1) hand3(1)];
        S=6;
        H=0.3;
        V=0.05;
end
        
% ----------------- %
% experimental data %
% ----------------- %
i=1;
for R=experimental_R
    file=['wall_profile' num2str(R) '.csv'];
    M=csvread(['./exp_data/' file],1);
    Z=M(1:S:end,1);
    U=M(1:S:end,2);
    plot(U(U>0),Z(U>0),['o' experimental_col(i)],'markersize',4);    
    [a,b]=max(U);
    text(1000*a,Z(b)+0.006,num2str(R),'color',experimental_col(i),'fontweight','bold');
    i=i+1;
end

% ---------------------- %
% simulation - developed %
% ---------------------- %
i=1;
for R=sim_dev_R 
  file=['r' num2str(R) 'D_U.csv'];
  M=csvread(['./sim_wall_profiles/Re375/developed/' file],1);
  Z=M(:,1)/0.001;
  U=M(:,3)/0.255; % Re375 - Um=0.255; Re800 - Um=0.544
  plot(U,Z,sim_dev_col(i),'linewidth',1);
  i=i+1;
end
i=1;
for R=-sim_dev_R 
  file=['r' num2str(abs(R)) 'D_U.csv'];
  M=csvread(['./sim_wall_profiles/Re800/developed/' file],1);
  Z=M(:,1)/0.001;
  U=M(:,3)/0.544;
  plot(-U,Z,sim_dev_col(i),'linewidth',1);
  i=i+1;
end

% ------------------- %
% Glauert - developed %
% ------------------- %
i=1;
M=csvread('Glauert profile.csv');
eta=M(:,1);
f=M(:,2);
for R=Glauert_R
  Z=8*eta*R^(5/4)/3^(3/4)/5^(1/4)/Re^(3/4);
  U=sqrt(15)*sqrt(Re)/32/R^(3/2)*2*f;
  plot(U, Z,'--','color',Glauert_col(i),'linewidth',1);
  i=i+1;
end

% ---------------------- %
% Schlichting - develped %
% ---------------------- %
i=1;
Z=0:0.002:0.1;
A=4.0;
for R=Schlichting_R
    U=A*R.*tanh(1.2*sqrt(A*Re)*Z);
    plot(U, Z,['--' Schlichting_col(i)],'linewidth',1);
    [a,b]=max(U);
    text(1000*a,Z(b),[' ' num2str(R)],'color',Schlichting_col(i),'fontweight','bold');
    i=i+1;
end


% --------------------- %
% Schlichting - uniform %
% --------------------- %
i=1;
Z=0:0.002:0.1;
A=1;
for R=-Schlichting_R
    U=-A*R.*tanh(1.2*sqrt(A*Re)*Z);
    plot(-U, Z,['--' Schlichting_col(i)],'linewidth',1);
    [a,b]=max(U);
    text(1000*a,Z(b),[' ' num2str(R)],'color',Schlichting_col(i),'fontweight','bold');
    i=i+1;
end


% ----------- %
% annotations %
% ----------- %
xlabel('\itu\rm [-]','fontangle','italic')
ylabel('\ity\rm [-]','fontangle','italic')
box on
legend(hands,'numerical','experimental','Glauert')
%text(H,V,'--- simulation')
%text(H,V-0.02,'- - Glauert')
%text(H,V-0.04,' o  experimental data')
for i=1:length(sim_dev_R)
   text(H,V-i*0.005,[num2str(sim_dev_R(i))],'color',sim_dev_col(i))
end
    
    
    