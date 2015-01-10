
clc;clear all;tic;
on=1;off=0;

FID=fopen('MDM_Output.dat', 'r');
A=fread(FID,4,'int32');

if exist('StoreM'),
    disp('Using current variables...');    
else
    X=A(1);
    Y=A(2);
    Z=A(3);
    Nr=A(4);
    
    ArrayMemSize=X*Y*Z;

    StoreM=zeros(X,Y,'double');

    disp('Reading main data...'); 
    
    popM=reshape(fread(FID,X*Y,'float32'),X,Y);

end

% Get Screen dimensions and set Main Window Dimensions
x = get(0,'ScreenSize'); ScreenDim=x(3:4);
MainWindowDim=[640 480];

% The graph window is initiated, with specified dimensions.
Figure1=figure('Position',[(ScreenDim-MainWindowDim)/2 MainWindowDim]);

Image=cat(3,(popM==0),(popM==2),zeros(X,Y));

subplot('position',[0.07 0.075 0.9 0.85]);

F1=imagesc(Image,[0 8]);
title('Green = Mussels; Red = Disturbed; Black = Bare','FontSize',12);  
colorbar; 
colormap('default'); axis image;drawnow;
    
for x=2:Z,
    popM=reshape(fread(FID,X*Y,'float32'),X,Y);
    Image=cat(3,(popM==0),(popM==2),zeros(X,Y));
    set(F1,'CData',Image);drawnow;
    
    set(Figure1,'Name',['Timestep ' num2str(x/Z*Nr) ' of ' num2str(Nr)]); 
    
end

fclose(FID);

toc

beep
