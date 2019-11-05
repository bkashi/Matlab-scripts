clear all; close all; clc;

%% load
global Nu Re Sd Hd Pr Cd
workpath='C:\Users\LENOVO\Desktop\Chapter 3 - Quantitative comparison of existing correlations\';
load([workpath 'experimental_and_correlated_data']);

%% reduce
minRe=158;   maxRe=100000;
minSd=2;     maxSd=31.1; 
minHd=1.6;   maxHd=20;
mintd=0.9;   maxtd=99;

i=1;
while i<=length(Nu)
            if Hd(i)<minHd || Hd(i)>maxHd || Sd(i)<minSd || Sd(i)>maxSd...
                 || Re(i)<minRe || Re(i)>maxRe || td(i)<mintd || td(i)>maxtd...
                 %|| distReturn(i)               
   
                 Re(i)=[]; Sd(i)=[]; Hd(i)=[]; Nu(i)=[]; Pr(i)=[];
                 d(i)=[];  N(i)=[]; td(i)=[]; Ar(i)=[]; reference(i)=[];
                 Figure(i)=[]; series(i)=[]; Meola(i)=[]; Robinson(i)=[];
                 Womac(i)=[]; Martin(i)=[];Cd(i)=[];
                 crossflow(i)=[]; uniformity(i)=[]; distReturn(i)=[];
                 L_(i)=[]; Re_L(i)=[]; Ar_(i)=[]; Lc(i)=[]; R(i)=[];
            else
                i=i+1;
            end
end


%%   optimize
UB=[100     0.5     0.9    1   1    1    1];
LB=[0.01    0.25    0.4   -1  -1   -1   -1];
options = optimoptions('ga','Display','diagnose','MaxGenerations',3000);
[X,fval] = ga(@RobinsonFitFcn,7,[],[],[],[],LB,UB,[],options);


%% plot
Nu_pred = X(1) * Pr.^X(2) .* Re.^X(3) .* Sd.^(X(4)+X(5)*(Hd>4)) .* Hd.^(X(6)+X(7)*(Hd>4));

loglog(Nu,Nu_pred,'.b'); hold on
set(gca,'Color',[1 1 1])
set(gca,'XTick',[1 10 100],'XTickLabel',[1 10 100],'YTick',[1 10 100],'YTickLabel',[1 10 100]);
plot([1 400],[1 400])

errThresh=0.25;
plot([0.01 900],[0.01 900]*(1+errThresh),':k')
plot([0.01 900],[0.01 900]*(1-errThresh),':k')

maxAx=max(Nu)*1.1;
minAx=4;
axis([minAx maxAx minAx maxAx])

text(maxAx/6,minAx*2,sprintf('%1.0f<Re<%1.0f\n %1.1f<S/d<%1.0f\n %1.1f<H/d<%1.0f\n %1.0f<t/d<%1.0f\n',...
                                    minRe,maxRe,minSd,maxSd,minHd,maxHd,mintd,maxtd),'Fontsize',8)
xlabel('Reported Nu_d','Fontsize',9)
ylabel('Calculated Nu_d','Fontsize',9)
title('Robinson, 2007','Fontsize',8)
hold on;
set(gca,'XScale','log','YScale','log','Box','on','Fontsize',8);
set(gca,'XTick',[1 10 100],'XTickLabel',[1 10 100],'YTick',[1 10 100],'YTickLabel',[1 10 100]);
plot([0.01 900],[0.01 900],'k')
plot([0.01 900],[0.01 900]*(1+errThresh),':k')
plot([0.01 900],[0.01 900]*(1-errThresh),':k')
axis([minAx maxAx minAx maxAx])

err=(Nu-Nu_pred)./Nu;
Naccurate=sum(abs(err)<errThresh,1);
mdl=fitlm(Nu,Nu_pred);
Rsq=mdl.Rsquared.Ordinary;
txt={['N=' num2str(length(Nu))] ;...
     ['N_{\Delta<25%}=' num2str(Naccurate/length(Nu)*100,'%1.0f') '%'] ;...
     ['\Delta_{avg.}=' num2str(mean(abs(err))*100,'%1.0f') '\pm' num2str(std(err)*100,'%1.0f') '%'] ;...
     ['SSE=' num2str(sum(err.^2),'%1.0f')] ;...
     ['SSE*=' num2str(sum(err.^2)/length(Nu),'%1.0f')] ;...
     ['r^2=' num2str(Rsq,2)]};
text(4.4,maxAx*0.34,txt,'Fontsize',8)
text(4.1,8.2,sprintf('+%1.0f%%',errThresh*100),'Fontsize',8)
text(6.8,4.7,sprintf('-%1.0f%%',errThresh*100),'Fontsize',8)
set(gcf,'Position',[180 60 340 290],'color','w')
%%
disp(num2str([[1:7]' X']))
