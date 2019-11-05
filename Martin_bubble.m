clear all
close all

r=0.25;
thetas=linspace(1/12*pi,11/12*pi,100);
Vs=pi*r^3./(sin(thetas)).^3.*(1-cos(thetas)).^2.*(2+cos(thetas));
plot(Vs,thetas/pi*180,'.')
xlabel('V')
ylabel('\theta (deg)')
close

figure


for V=0.01:0.01:1;
  theta=interp1(Vs,thetas,V);
  theta_deg=theta/pi*180;
  R=r/sin(theta);
  x=[];
  y=[];
  for phi=0:0.01:2*pi
    x(end+1)=R*cos(phi);
    y(end+1)=R*sin(phi)+R*sin(theta-pi/2);
  endfor
  plot(x(y>0),y(y>0),'-k');
  hold on
  plot([-2 2],[0 0],'-k');
  hold off
  ylim([0 1])
  xlim([-0.5 0.5])
  axis off
  axis square
  pause(0.1)
endfor

