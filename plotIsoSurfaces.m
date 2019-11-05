%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code reads smapled iso-surfaces in raw format,
% then sorts the point-pairs in ascending x-coordinate order,
% then reorders the points by proximinty, starting from x=0.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all previous data
clear all; close all;

% intialize plot and handels
figure; hold on
h1=plot(0,0);
h2=plot(0,0);
xlim([-0.0025 0.0025])
ylim([0 0.01])
box on;
set(gca,'DataAspectRatio',[1 1 1]);


% loop through all time folders
for i=0.0001:0.002:0.2000
  
  % convert time value to string
  folder=sprintf('%d',i);

  % read raw data file, skip header
    
    % the static contact angle case
    M1=dlmread(['stat/postProcessing/surfaces/' folder '/alpha.water_interface.raw'],' ',2,0);

    % the dynamic contact angle case
    M2=dlmread(['dyn/postProcessing/surfaces/' folder '/alpha.water_interface.raw'],' ',2,0);
    
  % remove z-column
  M1=M1(:,1:2);
  M2=M2(:,1:2);
  
  % sort by order
  M1=sortrows(M1);
  M2=sortrows(M2);
  
  % remove zeros
  M1(1:2,:)=[];
  M2(1:2,:)=[];
 
  % initialize N, the distance-ordered matrix
  N1=[];   N1(1,:)=M1(1,:);   j1=1;
  N2=[];   N2(1,:)=M2(1,:);   j2=1;
  
  % loop to reorder points of the static case
  while length(M1(:,1))>1

    % remove current point from original matrix
    M1(j1,:)=[];
    
    % calculate distances to current point
    distances1=M1-N1(end,:);
    
    % find index of closest point
    [~,j1]=min(sum(transpose(distances1.*distances1))');
    
    % add closest point to new matrix
    N1(end+1,:)=M1(j1,:);

  % end reordering loop
  endwhile
 
  % loop to reorder points of the dynamic case
  while length(M2(:,1))>1

    % remove current point from original matrix
    M2(j2,:)=[];
    
    % calculate distances to current point
    distances2=M2-N2(end,:);
    
    % find index of closest point
    [~,j2]=min(sum(transpose(distances2.*distances2))');
    
    % add closest point to new matrix
    N2(end+1,:)=M2(j2,:);

  % end reordering loop
  endwhile


 % update plot title
  title(['\itt\rm = ' folder 's']);
  
  % delete previous curves
  delete(h1);
  delete(h2);
  
  %plot x-y flipped points
  h1=plot([-N1(:,2); flipud(N1(:,2))],[N1(:,1); flipud(N1(:,1))],'-r','linewidth',2);
  h2=plot([-N2(:,2); flipud(N2(:,2))],[N2(:,1); flipud(N2(:,1))],'-g','linewidth',2);
  
  % update plot
  pause(0.01)
  
% end times loop  
endfor

 