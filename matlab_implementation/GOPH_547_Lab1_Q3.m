%GOPH 547 - Lab 1 - Forward Modelling for Gravity Anomalies 
%Safian Omar Qureshi
%ID: 10086638
%TA: Ye Sun
%Worked with Tian Yu, Ryan Younghwan Ok and TA 

%% Question 3:Part 1 - grid spacing=5m
clear; clc;

xm=[0,0,-10];        %mass anamoly centroid position
m=1.0*10^7;          %mass constant of anamoly
G=6.674*10^(-11);    %grav constant 

%set up meshgrid for a survey (5m spacing)
dx=5; xmin=-100; xmax=100;
dy=5; ymin=-100; ymax=100;

x=xmin:dx:xmax;
y=ymin:dy:ymax;

[X,Y]=meshgrid(x,y);
%%

%loop over all survey points for z=0 
Z0=0;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U0=zeros(size(X)); %potential function
gz0=zeros(size(X));%effect function

for j=1:nx
    for i=1:ny
        U0(i,j)=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm,m,G); %using functions 
        gz0(i,j)=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm,m,G);
    end 
end

%Loop over all survey points for z=10 (essentially
%identical algorithm to the loop defined above)
Z10=10;
ny=size(X,1); 
nx=size(X,2); 
U10=zeros(size(X));
gz10=zeros(size(X));

for j=1:nx
    for i=1:ny
        U10(i,j)=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm,m,G);
        gz10(i,j)=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm,m,G);
    end
end

%loop over all survey points for the z=100  
Z100=100;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U100=zeros(size(X));
gz100=zeros(size(X));

for j=1:nx
    for i=1:ny
        U100(i,j)=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm,m,G);
        gz100(i,j)=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm,m,G);
    end
end

%% Get min and max for U and gz colorbar limits

Umin=min(min(U100));
Umin=min(Umin,min(min(U0)));

Umax=max(max(U100));
Umax=max(Umax,max(max(U0)));


gzmin=min(min(gz100));
gzmin=min(gzmin,min(min(gz0)));

gzmax=max(max(gz100));
gzmax=max(gzmax,max(max(gz0)));

%% Creating figures for z=0,10 and 100 for meshgrid spacing 5

figure
subplot(3,2,1);
contourf(X,Y,U0);
hold on; % overlay grid points

plot(X,Y, 'xk', 'MarkerSize',2); %creating plot
axis equal; % plot x and y on equal scales
axis([xmin,xmax,ymin,ymax]); %define axis limits

xlabel('x [m]', 'fontweight','bold'); %x, y labels and title
ylabel('y [m]', 'fontweight','bold'); 
title('Gravitational Potential - depth 0m - 5m spacing')

h_c=colorbar; %add colorbar, change colormap at end 
ylabel(h_c, '[kg/s^2]','fontweight','bold');
caxis([Umin,Umax])
%%
subplot(3,2,2);
contourf(X,Y,gz0);
hold on; 

plot(X,Y, 'xk','MarkerSize',2);
axis equal; 
axis([xmin,xmax,ymin,ymax]); 

xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold'); 
title('Gravity Effect - depth 0m - 5m spacing')

h_c=colorbar;
ylabel(h_c, '[m/s^2]','fontweight','bold');
caxis([gzmin,gzmax])
%%
subplot(3,2,3)
contourf(X,Y,U10);
hold on; 

plot(X,Y, 'xk','MarkerSize',2);
axis equal;
axis([xmin,xmax,ymin,ymax]);
xlabel('x [m]','fontweight','bold'); 
ylabel('y [m]','fontweight','bold'); 
title('Gravitational Potential - depth 10m - 5m spacing')

h_c=colorbar; 
ylabel(h_c, '[kg/s^2]','fontweight','bold');
caxis([Umin,Umax])
%%
subplot(3,2,4);
contourf(X,Y,gz10);
hold on; 
plot(X,Y, 'xk','MarkerSize',2);
axis equal; 
axis([xmin,xmax,ymin,ymax]); 
xlabel('x [m]','fontweight','bold'); 
ylabel('y [m]','fontweight','bold'); 
title('Gravity Effect - depth 10m - 5m spacing')
h_c=colorbar; 
ylabel(h_c, '[m/s^2]','fontweight','bold');
caxis([gzmin,gzmax])

subplot(3,2,5)
contourf(X,Y,U100);
hold on; 
plot(X,Y, 'xk','MarkerSize',2);
axis equal; 
axis([xmin,xmax,ymin,ymax]); 
xlabel('x [m]','fontweight','bold'); 
ylabel('y [m]','fontweight','bold'); 
title('Gravitational Potential - depth 100m - 5m spacing')
h_c=colorbar; 
ylabel(h_c, '[kg/s^2]','fontweight','bold');
caxis([Umin,Umax])

subplot(3,2,6);
contourf(X,Y,gz100);
hold on; 
plot(X,Y, 'xk','MarkerSize',2);
axis equal; 
axis([xmin,xmax,ymin,ymax]); 
xlabel('x [m]','fontweight','bold'); 
ylabel('y [m]','fontweight','bold'); 
title('Gravity Effect - depth of 100m - 5m spacing')
h_c=colorbar; 
ylabel(h_c, '[m/s^2]','fontweight','bold');
caxis([gzmin,gzmax])

%changing color map 
cmap=colormap('parula');
cmap=flipud(cmap);
colormap(cmap);


%% Question 3:Part 2 - grid spacing=25m
clear; clc;

%nearly identical to script above, just changing gridspacing 

xm=[0,0,-10]; 
m=1.0*10^7; 
G=6.674*10^(-11);

%set up meshgrid for a survey (25m spacing)
dx=25; xmin=-100; xmax=100;
dy=25; ymin=-100; ymax=100;

x=xmin:dx:xmax;
y=ymin:dy:ymax;

[X,Y]=meshgrid(x,y);

%%

%simply copy pasted above script 
%loop over all survey points for z=0 
Z0=0;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U0=zeros(size(X)); %potential function
gz0=zeros(size(X));%effect function

for j=1:nx
    for i=1:ny
        U0(i,j)=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm,m,G); %using functions 
        gz0(i,j)=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm,m,G);
    end 
end

%Loop over all survey points for z=10 (essentially
%identical algorithm to the loop defined above)
Z10=10;
ny=size(X,1); 
nx=size(X,2); 
U10=zeros(size(X));
gz10=zeros(size(X));

for j=1:nx
    for i=1:ny
        U10(i,j)=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm,m,G);
        gz10(i,j)=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm,m,G);
    end
end

%loop over all survey points for the z=100  
Z100=100;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U100=zeros(size(X));
gz100=zeros(size(X));

for j=1:nx
    for i=1:ny
        U100(i,j)=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm,m,G);
        gz100(i,j)=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm,m,G);
    end
end

%% Get min and max for U and gz colorbar limits

Umin=min(min(U100));
Umin=min(Umin,min(min(U0)));

Umax=max(max(U100));
Umax=max(Umax,max(max(U0)));


gzmin=min(min(gz100));
gzmin=min(gzmin,min(min(gz0)));

gzmax=max(max(gz100));
gzmax=max(gzmax,max(max(gz0)));

%% Creating figures for z=0,10 and 100 for meshgrid spacing 25

figure
subplot(3,2,1);
contourf(X,Y,U0);
hold on; % overlay grid points

plot(X,Y, 'xk', 'MarkerSize',2); %creating plot
axis equal; % plot x and y on equal scales
axis([xmin,xmax,ymin,ymax]); %define axis limits

xlabel('x [m]', 'fontweight','bold'); %x, y labels and title
ylabel('y [m]', 'fontweight','bold'); 
title('Gravitational Potential - depth 0m - 25m spacing')

h_c=colorbar; %add colorbar, change colormap at end 
ylabel(h_c, '[kg/s^2]','fontweight','bold');
caxis([Umin,Umax])
%%
subplot(3,2,2);
contourf(X,Y,gz0);
hold on; 

plot(X,Y, 'xk','MarkerSize',2);
axis equal; 
axis([xmin,xmax,ymin,ymax]); 

xlabel('x [m]', 'fontweight','bold'); 
ylabel('y [m]', 'fontweight','bold'); 
title('Gravity Effect - depth 0m - 25m spacing')

h_c=colorbar;
ylabel(h_c, '[m/s^2]','fontweight','bold');
caxis([gzmin,gzmax])
%%
subplot(3,2,3)
contourf(X,Y,U10);
hold on; 

plot(X,Y, 'xk','MarkerSize',2);
axis equal;
axis([xmin,xmax,ymin,ymax]);
xlabel('x [m]','fontweight','bold'); 
ylabel('y [m]','fontweight','bold'); 
title('Gravitational Potential - depth 10m - 25m spacing')

h_c=colorbar; 
ylabel(h_c, '[kg/s^2]','fontweight','bold');
caxis([Umin,Umax])
%%
subplot(3,2,4);
contourf(X,Y,gz10);
hold on; 
plot(X,Y, 'xk','MarkerSize',2);
axis equal; 
axis([xmin,xmax,ymin,ymax]); 
xlabel('x [m]','fontweight','bold'); 
ylabel('y [m]','fontweight','bold'); 
title('Gravity Effect - depth 10m - 25m spacing')
h_c=colorbar; 
ylabel(h_c, '[m/s^2]','fontweight','bold');
caxis([gzmin,gzmax])

subplot(3,2,5)
contourf(X,Y,U100);
hold on; 
plot(X,Y, 'xk','MarkerSize',2);
axis equal; 
axis([xmin,xmax,ymin,ymax]); 
xlabel('x [m]','fontweight','bold'); 
ylabel('y [m]','fontweight','bold'); 
title('Gravitational Potential - depth 100m - 25m spacing')
h_c=colorbar; 
ylabel(h_c, '[kg/s^2]','fontweight','bold');
caxis([Umin,Umax])

subplot(3,2,6);
contourf(X,Y,gz100);
hold on; 
plot(X,Y, 'xk','MarkerSize',2);
axis equal; 
axis([xmin,xmax,ymin,ymax]); 
xlabel('x [m]','fontweight','bold'); 
ylabel('y [m]','fontweight','bold'); 
title('Gravity Effect - depth of 100m - 25m spacing')
h_c=colorbar; 
ylabel(h_c, '[m/s^2]','fontweight','bold');
caxis([gzmin,gzmax])

%changing color map 
cmap=colormap('parula');
cmap=flipud(cmap);
colormap(cmap);
