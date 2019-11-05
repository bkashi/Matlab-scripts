clear all; close all;

drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0, varargin{:} );

%-----------------------%
%  problem parameters   %
%-----------------------%

% Reynolds number
Re=1000;

% nozzle length (dimensionless)
l=100;

% nozzle-plate separation (dimensionless)
h=18;

% normalize nozzle-length & separation
L=l/Re;
H=h/Re;

%-----------------------%
%  draw nozzle & plate  %
%-----------------------%
if 1==1
  figure('position', [140 20 500 600]); hold on; box on; axis off;
  xlim([-6 6]); ylim([-1 21]);

  % paint fluid
  patch([-6 6 6 -6],[0 0 12 12],[1 1 1]*0.9,'edgecolor','none')
  patch([-1.2 1.2 1.2 -1.2],[12 12 18 18],[1 1 1]*0.9,'edgecolor','none')
  patch([-6 6 6 -6],[18 18 22 22],[1 1 1]*0.9,'edgecolor','none')
  
  % axis of symmetry
  plot([0 0],[-1 25],'-.k', 'linewidth',1.5)

  % plate
  plot([-6 6],[0 0],'-k','linewidth',4)

  % nozzle
  plot([-6 -1.2],[12 12],'-k','linewidth',2)
  plot([6 1.2],[12 12],'-k','linewidth',2)
  plot([-6 -1.2],[18 18],'-k','linewidth',2)
  plot([6 1.2],[12 12],'-k','linewidth',2)
  plot([6 1.2],[18 18],'-k','linewidth',2)
  plot([-1.2 -1.2],[12 18],'-k','linewidth',2)
  plot([1.2 1.2],[12 18],'-k','linewidth',2)
end

% annotations
text(1.4,1.4,'arriving profile','fontsize',12)
text(1.4,11,'issuinging profile','fontsize',12)
drawArrow([-4 -4],[6 12],'MaxHeadSize',0.15,'color','k')
drawArrow([-4 -4],[6 0],'MaxHeadSize',0.15,'color','k')
drawArrow([-4 -4],[15 18],'MaxHeadSize',0.3,'color','k')
drawArrow([-4 -4],[15 12],'MaxHeadSize',0.3,'color','k')
text(-4.3,14,['\itl\rm/\itd\rm=' num2str(l)],'fontname','times','fontsize',13,'rotation',90)
text(-4.3,5.5,['\ith\rm/\itd\rm=' num2str(h)],'fontname','times','fontsize',13,'rotation',90)
for x=-5.5:0.5:5.5
  drawArrow([x x],[-0.9 -0.15],'MaxHeadSize',1.2,'color','r','linewidth',2)
end
text(1,-1.2,'q','fontname','times','fontsize',20,'fontangle','italic','color','r')



%-------------------%
%  issuing profile  %
%-------------------%

y=linspace(0,1,100);
X=L*2;
delta=3.9611*(  2-exp(-750*X)/7-6*exp(-39*X)/7 -sqrt(2-exp(-750*X)/7-6*exp(-39*X)/7)  )  /  (2-exp(-750*X)/7-6*exp(-39*X)/7);
Uc=2-exp(-750*X)/7-6*exp(-39*X)/7;
if L<0.036
  lambda=2.99*(1-exp(-750*X)/7-6*exp(-39*X)/7)^0.39;
else
  lambda=0;
  delta=1;
end
U= (  2*y-2*y.^3+y.^4  +  lambda/6*(y-3*y.^2+3*y.^3-y.^4)  ) * Uc;
r=1-y*delta;

% add potential core
r1=0:0.01:r(end);
r=[r1 r(end:-1:1)];
U=[ones(1,length(r1))*U(end) U(end:-1:1)];

% jet cut-off
  [~,cutOff]=max(r(U>0.015*U(1)));
  
% plot issuing profile
if 1==1
  plot(r*1.2,12-U,'-b','linewidth',1)
  plot(-r*1.2,12-U,'-b','linewidth',1)
  plot([-1.2 1.2],[12 12],'-b','linewidth',0.8)
end

for q=1:24:cutOff
      drawArrow(1.2*[r(q) r(q)],[12  12-U(q)],'MaxHeadSize',0.6,'color','b','linewidth',2)
      drawArrow(-1.2*[r(q) r(q)],[12  12-U(q)],'MaxHeadSize',0.6,'color','b','linewidth',2)
end


%--------------------%
%       flight       %
%--------------------%

% rescale from r/R to r/d
r=r/2;

% domain size
b=2.0;

% number or required roots
N=64;

% find roots - beta values
t=linspace(0,10*N,10000);
B=besselj(0,t*b);
j=1;
for i=1:N
  while sign(B(j))*sign(B(j+1))==1
    j=j+1;
  end
  beta(i)=(t(j)+t(j+1))/2;
  j=j+1;
end

% integrate velocity profile
for i=1:N
  I(i)=trapz(r,r.*besselj(0,beta(i)*r).*U);
end

% radial domain
r=linspace(0,b);

% flight distance
X=[1/4 2/4 3/4 4/4]*H/U(1);

for n=1:length(X)
  
  % initialize u
  u=zeros(1,length(r));

  % sum u terms
  for i=1:N
    u=u+ 2/b^2 * exp(-(beta(i))^2*X(n)) * besselj(0,beta(i)*r) / (besselj(1,beta(i)*b))^2 * I(i); 
  end
    
  % jet cut-off
  [~,cutOff]=max(r(u>0.015*u(1)));
  
  % plot u up to jet width
  plot(2.4*r(1:cutOff),9.8-(n-1)*10/4-u(1:cutOff),'-b','linewidth',1) 
  plot(-2.4*r(1:cutOff),9.8-(n-1)*10/4-u(1:cutOff),'-b','linewidth',1)
  plot([-2.4*r(cutOff) 2.4*r(cutOff)],[9.8-(n-1)*10/4 9.8-(n-1)*10/4],'-b','linewidth',0.8)
  
  for q=1:6:cutOff
      drawArrow(2.4*[r(q) r(q)],[9.8-(n-1)*10/4  9.8-(n-1)*10/4-u(q)],'MaxHeadSize',0.6,'color','b','linewidth',2)
      drawArrow(-2.4*[r(q) r(q)],[9.8-(n-1)*10/4  9.8-(n-1)*10/4-u(q)],'MaxHeadSize',0.6,'color','b','linewidth',2)
  end
  
end  


%-------------------------%
% local velocity gradient %
%-------------------------%

% radial velocity gradient
Rw=min( r(u<=0.81*u(1)) );
A0=0.88*u(1)/2/Rw;

% stagnation zone height
Zw=(A0/0.88/U(1))^(0.77/2)/U(1)^0.77;




