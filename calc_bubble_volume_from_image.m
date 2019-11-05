clear all; close all;

% prepare sums collection for all images
I=[];

% loop through all image files
for f=[2140]

    % convert filename to string
    F=sprintf('%1.6f.png',f/1e6);
    
    % read image file to matrix
    M=imread(F);
    
    % invert image
    M=255-M;
    
    % stretch grayscale to B&W
    M(M<128)=0;
    M(M>128)=255;
    
    % find axis or symmetry
    Y0=[];
    for x=1:size(M,2)
        y0=max(find(M(:,x,1)==255));
        if size(y0)>0 Y0(end+1)=y0; end
    end
    
    % mark axis in red on dimensions 2 & 3
    M(max(Y0),:,2:3)=0;
    
    % trim image below axis
    M=M(1:max(Y0),:,:);
    
    % flip image such that the axis lays on y=0
    M(:,:,1)=flipud(M(:,:,1));
    M(:,:,2)=flipud(M(:,:,2));
    M(:,:,3)=flipud(M(:,:,3));
    
    % prepare a radial distance matrix for volume integration
    R=[1:length(M(:,1,1))]';
    R=repmat(R,1,length(M(1,:,1)));
    
    % reduce M to 2D
    %M=double(M(:,:,1))/255;
    
    % mutipliy matrices
    Q=(double(M(:,:,1))/255).*R;
    
    % sum Q
    I(end+1)=round(2*pi*sum(sum(Q)));      

    % mean equivalent sphere
    r=int16(mean((3/4/pi*I).^(1/3)));

    % draw equivalent circle on image dimension 1

        % find circle center in x-direction
        x0=int16(length(M(1,:,1))/2);
        
        % define x vector
        x=int16((x0-r):(x0+r));
        
        % caqlculate y vector
        y=int16(sqrt(r^2-(x-x0).^2));
        
        % "plot" circle into image
        for i=1:length(x)
            M(y(i)+1,x(i),1)=255;
            M(y(i)+1,x(i),2:3)=0;
        end
        
    % output ratio of peak point to r
    r/maxmax(M)
    
    % show & save image
    figure
    imshow(M); 
    imwrite(M,[F(1:end-4) '_.png'])    

end

% min. to max. volume ratio
max(I)/min(I) 
    
    