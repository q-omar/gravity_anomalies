%GOPH 547 - Lab 1 - Forward Modelling for Gravity Anomalies 
%Safian Omar Qureshi
%ID: 10086638
%TA: Ye Sun
%Worked with Tian Yu, Ryan Younghwan Ok and TA 

%% Question 4a and 4b  
clear; clc;

m=1.0*10^7; %Total Mass of 5 anomalies

%normal distribution constraints

mu_m=m/5;
mu_x=0;
mu_y=0;
mu_z=-10;

sig_m=m/100;
sig_x=20;
sig_y=20;
sig_z=2;

%% Generate 4 random values for set 1
m_rand=randn(4,1)*sig_m+mu_m;
x_rand=randn(4,1)*sig_x+mu_x;
y_rand=randn(4,1)*sig_y+mu_y;
z_rand=randn(4,1)*sig_z+mu_z;

% Calculate the 5th value for set 1
m_5=m-sum(m_rand);
x_5=(-sum(m_rand.*x_rand))/m_5;
y_5=(-sum(m_rand.*y_rand))/m_5;
z_5=(mu_z*m-(sum(m_rand.*z_rand)))/m_5;

save mass_set_1 %generating mass set 1 

%% Generate 4 random values for set 2
m_rand=randn(4,1)*sig_m+mu_m;
x_rand=randn(4,1)*sig_x+mu_x;
y_rand=randn(4,1)*sig_y+mu_y;
z_rand=randn(4,1)*sig_z+mu_z;

% Calculate the 5th value for set 2 
m_5=m-sum(m_rand);
x_5=(-sum(m_rand.*x_rand))/m_5;
y_5=(-sum(m_rand.*y_rand))/m_5;
z_5=(mu_z*m-(sum(m_rand.*z_rand)))/m_5;

save mass_set_2 %generating mass set 2 


%% Generate 4 random values for set 3
m_rand=randn(4,1)*sig_m+mu_m;
x_rand=randn(4,1)*sig_x+mu_x;
y_rand=randn(4,1)*sig_y+mu_y;
z_rand=randn(4,1)*sig_z+mu_z;

% Calculate the 5th value for set 3
m_5=m-sum(m_rand);
x_5=(-sum(m_rand.*x_rand))/m_5;
y_5=(-sum(m_rand.*y_rand))/m_5;
z_5=(mu_z*m-(sum(m_rand.*z_rand)))/m_5;

save mass_set_3 %generating mass set 3

%% Question 4c mass set 1 - grid spacing 5m

clear; clc;

load mass_set_1

G=6.674*10^(-11);

%Determining position of each individual random mass point
xm_1=[x_rand(1),y_rand(1),z_rand(1)];
xm_2=[x_rand(2),y_rand(2),z_rand(2)];
xm_3=[x_rand(3),y_rand(3),z_rand(3)];
xm_4=[x_rand(4),y_rand(4),z_rand(4)];
xm_5=[x_5,y_5,z_5];%Location of 5th mass is calculated previously

%Reading the mass of each anomaly  
m_1=m_rand(1);
m_2=m_rand(2);
m_3=m_rand(3);
m_4=m_rand(4);


%%
%Setting up survey meshgrid as in question 3: grid spacing 5m
dx=5; xmin=-100; xmax=100;
dy=5; ymin=-100; ymax=100;

x=xmin:dx:xmax;
y=ymin:dy:ymax;

[X,Y]=meshgrid(x,y);

%%
%loop over all survey points for z=0 as question 3 
Z0=0;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U0=zeros(size(X));%potential function
gz0=zeros(size(X));%effect function 

for j=1:nx
    for i=1:ny
        U0_1=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_1,m_1,G);
        U0_2=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_2,m_2,G);
        U0_3=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_3,m_3,G);
        U0_4=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_4,m_4,G);
        U0_5=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_5,m_5,G);
        U0(i,j)=U0_1+U0_2+U0_3+U0_4+U0_5; %summing contributions from each of the individual masses 
        gz0_1=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_1,m_1,G);
        gz0_2=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_2,m_2,G);
        gz0_3=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_3,m_3,G);
        gz0_4=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_4,m_4,G);
        gz0_5=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_5,m_5,G);
        gz0(i,j)=gz0_1+gz0_2+gz0_3+gz0_4+gz0_5; %summing contributions from each of the individual masses
    end
end

%loop over all survey points for z=10 
Z10=10;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U10=zeros(size(X));
gz10=zeros(size(X));

for j=1:nx
    for i=1:ny
        U10_1=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_1,m_1,G);
        U10_2=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_2,m_2,G);
        U10_3=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_3,m_3,G);
        U10_4=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_4,m_4,G);
        U10_5=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_5,m_5,G);
        U10(i,j)=U10_1+U10_2+U10_3+U10_4+U10_5;
        gz10_1=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_1,m_1,G);
        gz10_2=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_2,m_2,G);
        gz10_3=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_3,m_3,G);
        gz10_4=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_4,m_4,G);
        gz10_5=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_5,m_5,G);
        gz10(i,j)=gz10_1+gz10_2+gz10_3+gz10_4+gz10_5;
    end
end

%loop over all survey points for z=10 
Z100=100;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U100=zeros(size(X));
gz100=zeros(size(X));

for j=1:nx
    for i=1:ny
        U100_1=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_1,m_1,G);
        U100_2=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_2,m_2,G);
        U100_3=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_3,m_3,G);
        U100_4=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_4,m_4,G);
        U100_5=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_5,m_5,G);
        U100(i,j)=U100_1+U100_2+U100_3+U100_4+U100_5;
        gz100_1=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_1,m_1,G);
        gz100_2=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_2,m_2,G);
        gz100_3=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_3,m_3,G);
        gz100_4=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_4,m_4,G);
        gz100_5=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_5,m_5,G);
        gz100(i,j)=gz100_1+gz100_2+gz100_3+gz100_4+gz100_5;
    end
end

%% get min and max for U and gz for colorbar limits
Umin=min(min(U100));
Umin=min(Umin,min(min(U0)));

Umax=max(max(U100));
Umax=max(Umax,max(max(U0)));

gzmin=min(min(gz100));
gzmin=min(gzmin,min(min(gz0)));

gzmax=max(max(gz100));
gzmax=max(gzmax,max(max(gz0)));
%% Creating figures contour plots for z=0,10 and 100 for meshgrid spacing 5

figure
subplot(3,2,1);
contourf(X,Y,U0);
hold on; % overlay grid points

plot(X,Y, 'xk', 'MarkerSize',2); %creating plot
axis equal; % plot x and y on equal scales
axis([xmin,xmax,ymin,ymax]); %define axis limits

xlabel('x [m]', 'fontweight','bold'); %x, y labels and title
ylabel('y [m]', 'fontweight','bold'); 
title('Gravitational Potential - depth 0m - 5m spacing - 5 random masses')

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
title('Gravitational Effect - depth 0m - 5m spacing - 5 random masses')

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
title('Gravitational Potential - depth 10m - 5m spacing - 5 random masses (set1)')

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
title('Gravitational Effect - depth 10m - 5m spacing - 5 random masses (set1)')
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
title('Gravitational Potential - depth 100m - 5m spacing - 5 random masses (set1)')
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
title('Gravitational Effect - depth 0m - 5m spacing - 5 random masses (set1)')
h_c=colorbar; 
ylabel(h_c, '[m/s^2]','fontweight','bold');
caxis([gzmin,gzmax])

cmap=colormap('parula');
cmap=flipud(cmap);
colormap(cmap);
%%

%Setting up survey meshgrid for 25m 
dx=25; xmin=-100; xmax=100;
dy=25; ymin=-100; ymax=100;

x=xmin:dx:xmax;
y=ymin:dy:ymax;

[X,Y]=meshgrid(x,y);

%Loop over survey pts for the z=0 case. 
%loop over all survey points for z=0 as question 3 
Z0=0;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U0=zeros(size(X));%potential function
gz0=zeros(size(X));%effect function 

for j=1:nx
    for i=1:ny
        U0_1=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_1,m_1,G);
        U0_2=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_2,m_2,G);
        U0_3=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_3,m_3,G);
        U0_4=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_4,m_4,G);
        U0_5=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_5,m_5,G);
        U0(i,j)=U0_1+U0_2+U0_3+U0_4+U0_5; %summing contributions from each of the individual masses 
        gz0_1=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_1,m_1,G);
        gz0_2=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_2,m_2,G);
        gz0_3=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_3,m_3,G);
        gz0_4=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_4,m_4,G);
        gz0_5=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_5,m_5,G);
        gz0(i,j)=gz0_1+gz0_2+gz0_3+gz0_4+gz0_5; %summing contributions from each of the individual masses
    end
end

%loop over all survey points for z=10 
Z10=10;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U10=zeros(size(X));
gz10=zeros(size(X));

for j=1:nx
    for i=1:ny
        U10_1=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_1,m_1,G);
        U10_2=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_2,m_2,G);
        U10_3=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_3,m_3,G);
        U10_4=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_4,m_4,G);
        U10_5=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_5,m_5,G);
        U10(i,j)=U10_1+U10_2+U10_3+U10_4+U10_5;
        gz10_1=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_1,m_1,G);
        gz10_2=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_2,m_2,G);
        gz10_3=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_3,m_3,G);
        gz10_4=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_4,m_4,G);
        gz10_5=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_5,m_5,G);
        gz10(i,j)=gz10_1+gz10_2+gz10_3+gz10_4+gz10_5;
    end
end

%loop over all survey points for z=100
Z100=100;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U100=zeros(size(X));
gz100=zeros(size(X));

for j=1:nx
    for i=1:ny
        U100_1=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_1,m_1,G);
        U100_2=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_2,m_2,G);
        U100_3=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_3,m_3,G);
        U100_4=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_4,m_4,G);
        U100_5=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_5,m_5,G);
        U100(i,j)=U100_1+U100_2+U100_3+U100_4+U100_5;
        gz100_1=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_1,m_1,G);
        gz100_2=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_2,m_2,G);
        gz100_3=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_3,m_3,G);
        gz100_4=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_4,m_4,G);
        gz100_5=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_5,m_5,G);
        gz100(i,j)=gz100_1+gz100_2+gz100_3+gz100_4+gz100_5;
    end
end

%% get min and max for U and gz for colorbar limits
Umin=min(min(U100));
Umin=min(Umin,min(min(U0)));

Umax=max(max(U100));
Umax=max(Umax,max(max(U0)));

gzmin=min(min(gz100));
gzmin=min(gzmin,min(min(gz0)));

gzmax=max(max(gz100));
gzmax=max(gzmax,max(max(gz0)));
%% Creating figures contour plots for z=0,10 and 100 for meshgrid spacing 25

figure
subplot(3,2,1);
contourf(X,Y,U0);
hold on; % overlay grid points

plot(X,Y, 'xk', 'MarkerSize',2); %creating plot
axis equal; % plot x and y on equal scales
axis([xmin,xmax,ymin,ymax]); %define axis limits

xlabel('x [m]', 'fontweight','bold'); %x, y labels and title
ylabel('y [m]', 'fontweight','bold'); 
title('Gravitational Potential - depth 0m - 25m spacing - 5 random masses (set1)')

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
title('Gravitational Effect - depth 0m - 25m spacing - 5 random masses (set1)')

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
title('Gravitational Potential - depth 10m - 25m spacing - 5 random masses (set1)')

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
title('Gravitational Effect - depth 10m - 25m spacing - 5 random masses (set1)')
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
title('Gravitational Potential - depth 100m - 25m spacing - 5 random masses (set1)')
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
title('Gravitational E - depth 0m - 25m spacing - 5 random masses (set1)')
h_c=colorbar; 
ylabel(h_c, '[m/s^2]','fontweight','bold');
caxis([gzmin,gzmax])

cmap=colormap('parula');
cmap=flipud(cmap);
colormap(cmap);

%%

%% Question 4c: Mass set 2 - grid spacing=5m
clear; clc;
G=6.674*10^(-11);
load mass_set_2

G=6.674*10^(-11);

%Determining position of each individual random mass point
xm_1=[x_rand(1),y_rand(1),z_rand(1)];
xm_2=[x_rand(2),y_rand(2),z_rand(2)];
xm_3=[x_rand(3),y_rand(3),z_rand(3)];
xm_4=[x_rand(4),y_rand(4),z_rand(4)];
xm_5=[x_5,y_5,z_5];%Location of 5th mass is calculated previously

%Reading the mass of each anomaly  
m_1=m_rand(1);
m_2=m_rand(2);
m_3=m_rand(3);
m_4=m_rand(4);


%%
%Setting up survey meshgrid as in question 3: grid spacing 5m
dx=5; xmin=-100; xmax=100;
dy=5; ymin=-100; ymax=100;

x=xmin:dx:xmax;
y=ymin:dy:ymax;

[X,Y]=meshgrid(x,y);

%%
%loop over all survey points for z=0 as question 3 
Z0=0;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U0=zeros(size(X));%potential function
gz0=zeros(size(X));%effect function 

for j=1:nx
    for i=1:ny
        U0_1=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_1,m_1,G);
        U0_2=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_2,m_2,G);
        U0_3=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_3,m_3,G);
        U0_4=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_4,m_4,G);
        U0_5=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_5,m_5,G);
        U0(i,j)=U0_1+U0_2+U0_3+U0_4+U0_5; %summing contributions from each of the individual masses 
        gz0_1=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_1,m_1,G);
        gz0_2=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_2,m_2,G);
        gz0_3=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_3,m_3,G);
        gz0_4=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_4,m_4,G);
        gz0_5=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_5,m_5,G);
        gz0(i,j)=gz0_1+gz0_2+gz0_3+gz0_4+gz0_5; %summing contributions from each of the individual masses
    end
end

%loop over all survey points for z=10 
Z10=10;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U10=zeros(size(X));
gz10=zeros(size(X));

for j=1:nx
    for i=1:ny
        U10_1=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_1,m_1,G);
        U10_2=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_2,m_2,G);
        U10_3=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_3,m_3,G);
        U10_4=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_4,m_4,G);
        U10_5=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_5,m_5,G);
        U10(i,j)=U10_1+U10_2+U10_3+U10_4+U10_5;
        gz10_1=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_1,m_1,G);
        gz10_2=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_2,m_2,G);
        gz10_3=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_3,m_3,G);
        gz10_4=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_4,m_4,G);
        gz10_5=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_5,m_5,G);
        gz10(i,j)=gz10_1+gz10_2+gz10_3+gz10_4+gz10_5;
    end
end

%loop over all survey points for z=10 
Z100=100;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U100=zeros(size(X));
gz100=zeros(size(X));

for j=1:nx
    for i=1:ny
        U100_1=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_1,m_1,G);
        U100_2=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_2,m_2,G);
        U100_3=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_3,m_3,G);
        U100_4=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_4,m_4,G);
        U100_5=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_5,m_5,G);
        U100(i,j)=U100_1+U100_2+U100_3+U100_4+U100_5;
        gz100_1=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_1,m_1,G);
        gz100_2=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_2,m_2,G);
        gz100_3=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_3,m_3,G);
        gz100_4=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_4,m_4,G);
        gz100_5=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_5,m_5,G);
        gz100(i,j)=gz100_1+gz100_2+gz100_3+gz100_4+gz100_5;
    end
end

%% get min and max for U and gz for colorbar limits
Umin=min(min(U100));
Umin=min(Umin,min(min(U0)));

Umax=max(max(U100));
Umax=max(Umax,max(max(U0)));

gzmin=min(min(gz100));
gzmin=min(gzmin,min(min(gz0)));

gzmax=max(max(gz100));
gzmax=max(gzmax,max(max(gz0)));
%% Creating figures contour plots for z=0,10 and 100 for meshgrid spacing 5

figure
subplot(3,2,1);
contourf(X,Y,U0);
hold on; % overlay grid points

plot(X,Y, 'xk', 'MarkerSize',2); %creating plot
axis equal; % plot x and y on equal scales
axis([xmin,xmax,ymin,ymax]); %define axis limits

xlabel('x [m]', 'fontweight','bold'); %x, y labels and title
ylabel('y [m]', 'fontweight','bold'); 
title('Gravitational Potential - depth 0m - 5m spacing - 5 random masses (set2)')

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
title('Gravitational Effect - depth 0m - 5m spacing - 5 random masses (set2)')

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
title('Gravitational Potential - depth 10m - 5m spacing - 5 random masses (set2)')

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
title('Gravitational Effect - depth 10m - 5m spacing - 5 random masses (set2)')
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
title('Gravitational Potential - depth 100m - 5m spacing - 5 random masses (set2)')
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
title('Gravitational Effect - depth 0m - 5m spacing - 5 random masses (set2)')
h_c=colorbar; 
ylabel(h_c, '[m/s^2]','fontweight','bold');
caxis([gzmin,gzmax])

cmap=colormap('parula');
cmap=flipud(cmap);
colormap(cmap);
%%

%Setting up survey meshgrid for 25m 
dx=25; xmin=-100; xmax=100;
dy=25; ymin=-100; ymax=100;

x=xmin:dx:xmax;
y=ymin:dy:ymax;

[X,Y]=meshgrid(x,y);

%Loop over survey pts for the z=0 case. 
%loop over all survey points for z=0 as question 3 
Z0=0;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U0=zeros(size(X));%potential function
gz0=zeros(size(X));%effect function 

for j=1:nx
    for i=1:ny
        U0_1=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_1,m_1,G);
        U0_2=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_2,m_2,G);
        U0_3=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_3,m_3,G);
        U0_4=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_4,m_4,G);
        U0_5=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_5,m_5,G);
        U0(i,j)=U0_1+U0_2+U0_3+U0_4+U0_5; %summing contributions from each of the individual masses 
        gz0_1=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_1,m_1,G);
        gz0_2=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_2,m_2,G);
        gz0_3=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_3,m_3,G);
        gz0_4=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_4,m_4,G);
        gz0_5=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_5,m_5,G);
        gz0(i,j)=gz0_1+gz0_2+gz0_3+gz0_4+gz0_5; %summing contributions from each of the individual masses
    end
end

%loop over all survey points for z=10 
Z10=10;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U10=zeros(size(X));
gz10=zeros(size(X));

for j=1:nx
    for i=1:ny
        U10_1=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_1,m_1,G);
        U10_2=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_2,m_2,G);
        U10_3=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_3,m_3,G);
        U10_4=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_4,m_4,G);
        U10_5=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_5,m_5,G);
        U10(i,j)=U10_1+U10_2+U10_3+U10_4+U10_5;
        gz10_1=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_1,m_1,G);
        gz10_2=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_2,m_2,G);
        gz10_3=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_3,m_3,G);
        gz10_4=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_4,m_4,G);
        gz10_5=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_5,m_5,G);
        gz10(i,j)=gz10_1+gz10_2+gz10_3+gz10_4+gz10_5;
    end
end

%loop over all survey points for z=100
Z100=100;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U100=zeros(size(X));
gz100=zeros(size(X));

for j=1:nx
    for i=1:ny
        U100_1=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_1,m_1,G);
        U100_2=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_2,m_2,G);
        U100_3=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_3,m_3,G);
        U100_4=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_4,m_4,G);
        U100_5=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_5,m_5,G);
        U100(i,j)=U100_1+U100_2+U100_3+U100_4+U100_5;
        gz100_1=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_1,m_1,G);
        gz100_2=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_2,m_2,G);
        gz100_3=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_3,m_3,G);
        gz100_4=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_4,m_4,G);
        gz100_5=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_5,m_5,G);
        gz100(i,j)=gz100_1+gz100_2+gz100_3+gz100_4+gz100_5;
    end
end

%% get min and max for U and gz for colorbar limits
Umin=min(min(U100));
Umin=min(Umin,min(min(U0)));

Umax=max(max(U100));
Umax=max(Umax,max(max(U0)));

gzmin=min(min(gz100));
gzmin=min(gzmin,min(min(gz0)));

gzmax=max(max(gz100));
gzmax=max(gzmax,max(max(gz0)));
%% Creating figures contour plots for z=0,10 and 100 for meshgrid spacing 25

figure
subplot(3,2,1);
contourf(X,Y,U0);
hold on; % overlay grid points

plot(X,Y, 'xk', 'MarkerSize',2); %creating plot
axis equal; % plot x and y on equal scales
axis([xmin,xmax,ymin,ymax]); %define axis limits

xlabel('x [m]', 'fontweight','bold'); %x, y labels and title
ylabel('y [m]', 'fontweight','bold'); 
title('Gravitational Potential - depth 0m - 25m spacing - 5 random masses (set2)')

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
title('Gravitational Effect - depth 0m - 25m spacing - 5 random masses (set2)')

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
title('Gravitational Potential - depth 10m - 25m spacing - 5 random masses (set2)')

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
title('Gravitational Effect - depth 10m - 25m spacing - 5 random masses (set2)')
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
title('Gravitational Potential - depth 100m - 25m spacing - 5 random masses (set2)')
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
title('Gravitational E - depth 0m - 25m spacing - 5 random masses (set2)')
h_c=colorbar; 
ylabel(h_c, '[m/s^2]','fontweight','bold');
caxis([gzmin,gzmax])

cmap=colormap('parula');
cmap=flipud(cmap);
colormap(cmap);

%%
%% Question 4c: Mass set 3 - grid spacing=5m
clear; clc;

load mass_set_3

G=6.674*10^(-11);

%Determining position of each individual random mass point
xm_1=[x_rand(1),y_rand(1),z_rand(1)];
xm_2=[x_rand(2),y_rand(2),z_rand(2)];
xm_3=[x_rand(3),y_rand(3),z_rand(3)];
xm_4=[x_rand(4),y_rand(4),z_rand(4)];
xm_5=[x_5,y_5,z_5];%Location of 5th mass is calculated previously

%Reading the mass of each anomaly  
m_1=m_rand(1);
m_2=m_rand(2);
m_3=m_rand(3);
m_4=m_rand(4);


%%
%Setting up survey meshgrid as in question 3: grid spacing 5m
dx=5; xmin=-100; xmax=100;
dy=5; ymin=-100; ymax=100;

x=xmin:dx:xmax;
y=ymin:dy:ymax;

[X,Y]=meshgrid(x,y);

%%
%loop over all survey points for z=0 as question 3 
Z0=0;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U0=zeros(size(X));%potential function
gz0=zeros(size(X));%effect function 

for j=1:nx
    for i=1:ny
        U0_1=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_1,m_1,G);
        U0_2=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_2,m_2,G);
        U0_3=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_3,m_3,G);
        U0_4=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_4,m_4,G);
        U0_5=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_5,m_5,G);
        U0(i,j)=U0_1+U0_2+U0_3+U0_4+U0_5; %summing contributions from each of the individual masses 
        gz0_1=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_1,m_1,G);
        gz0_2=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_2,m_2,G);
        gz0_3=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_3,m_3,G);
        gz0_4=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_4,m_4,G);
        gz0_5=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_5,m_5,G);
        gz0(i,j)=gz0_1+gz0_2+gz0_3+gz0_4+gz0_5; %summing contributions from each of the individual masses
    end
end

%loop over all survey points for z=10 
Z10=10;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U10=zeros(size(X));
gz10=zeros(size(X));

for j=1:nx
    for i=1:ny
        U10_1=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_1,m_1,G);
        U10_2=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_2,m_2,G);
        U10_3=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_3,m_3,G);
        U10_4=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_4,m_4,G);
        U10_5=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_5,m_5,G);
        U10(i,j)=U10_1+U10_2+U10_3+U10_4+U10_5;
        gz10_1=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_1,m_1,G);
        gz10_2=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_2,m_2,G);
        gz10_3=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_3,m_3,G);
        gz10_4=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_4,m_4,G);
        gz10_5=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_5,m_5,G);
        gz10(i,j)=gz10_1+gz10_2+gz10_3+gz10_4+gz10_5;
    end
end

%loop over all survey points for z=10 
Z100=100;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U100=zeros(size(X));
gz100=zeros(size(X));

for j=1:nx
    for i=1:ny
        U100_1=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_1,m_1,G);
        U100_2=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_2,m_2,G);
        U100_3=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_3,m_3,G);
        U100_4=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_4,m_4,G);
        U100_5=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_5,m_5,G);
        U100(i,j)=U100_1+U100_2+U100_3+U100_4+U100_5;
        gz100_1=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_1,m_1,G);
        gz100_2=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_2,m_2,G);
        gz100_3=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_3,m_3,G);
        gz100_4=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_4,m_4,G);
        gz100_5=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_5,m_5,G);
        gz100(i,j)=gz100_1+gz100_2+gz100_3+gz100_4+gz100_5;
    end
end

%% get min and max for U and gz for colorbar limits
Umin=min(min(U100));
Umin=min(Umin,min(min(U0)));

Umax=max(max(U100));
Umax=max(Umax,max(max(U0)));

gzmin=min(min(gz100));
gzmin=min(gzmin,min(min(gz0)));

gzmax=max(max(gz100));
gzmax=max(gzmax,max(max(gz0)));
%% Creating figures contour plots for z=0,10 and 100 for meshgrid spacing 5

figure
subplot(3,2,1);
contourf(X,Y,U0);
hold on; % overlay grid points

plot(X,Y, 'xk', 'MarkerSize',2); %creating plot
axis equal; % plot x and y on equal scales
axis([xmin,xmax,ymin,ymax]); %define axis limits

xlabel('x [m]', 'fontweight','bold'); %x, y labels and title
ylabel('y [m]', 'fontweight','bold'); 
title('Gravitational Potential - depth 0m - 5m spacing - 5 random masses (set3)')

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
title('Gravitational Effect - depth 0m - 5m spacing - 5 random masses (set3)')

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
title('Gravitational Potential - depth 10m - 5m spacing - 5 random masses (set3)')

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
title('Gravitational Effect - depth 10m - 5m spacing - 5 random masses (set3)')
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
title('Gravitational Potential - depth 100m - 5m spacing - 5 random masses (set3)')
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
title('Gravitational Effect - depth 0m - 5m spacing - 5 random masses (set3)')
h_c=colorbar; 
ylabel(h_c, '[m/s^2]','fontweight','bold');
caxis([gzmin,gzmax])

cmap=colormap('parula');
cmap=flipud(cmap);
colormap(cmap);
%%

%Setting up survey meshgrid for 25m 
dx=25; xmin=-100; xmax=100;
dy=25; ymin=-100; ymax=100;

x=xmin:dx:xmax;
y=ymin:dy:ymax;

[X,Y]=meshgrid(x,y);

%Loop over survey pts for the z=0 case. 
%loop over all survey points for z=0 as question 3 
Z0=0;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U0=zeros(size(X));%potential function
gz0=zeros(size(X));%effect function 

for j=1:nx
    for i=1:ny
        U0_1=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_1,m_1,G);
        U0_2=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_2,m_2,G);
        U0_3=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_3,m_3,G);
        U0_4=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_4,m_4,G);
        U0_5=grav_pot_point(([X(i,j),Y(i,j),Z0]),xm_5,m_5,G);
        U0(i,j)=U0_1+U0_2+U0_3+U0_4+U0_5; %summing contributions from each of the individual masses 
        gz0_1=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_1,m_1,G);
        gz0_2=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_2,m_2,G);
        gz0_3=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_3,m_3,G);
        gz0_4=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_4,m_4,G);
        gz0_5=grav_eff_point(([X(i,j),Y(i,j),Z0]),xm_5,m_5,G);
        gz0(i,j)=gz0_1+gz0_2+gz0_3+gz0_4+gz0_5; %summing contributions from each of the individual masses
    end
end

%loop over all survey points for z=10 
Z10=10;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U10=zeros(size(X));
gz10=zeros(size(X));

for j=1:nx
    for i=1:ny
        U10_1=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_1,m_1,G);
        U10_2=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_2,m_2,G);
        U10_3=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_3,m_3,G);
        U10_4=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_4,m_4,G);
        U10_5=grav_pot_point(([X(i,j),Y(i,j),Z10]),xm_5,m_5,G);
        U10(i,j)=U10_1+U10_2+U10_3+U10_4+U10_5;
        gz10_1=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_1,m_1,G);
        gz10_2=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_2,m_2,G);
        gz10_3=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_3,m_3,G);
        gz10_4=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_4,m_4,G);
        gz10_5=grav_eff_point(([X(i,j),Y(i,j),Z10]),xm_5,m_5,G);
        gz10(i,j)=gz10_1+gz10_2+gz10_3+gz10_4+gz10_5;
    end
end

%loop over all survey points for z=100
Z100=100;
ny=size(X,1); % # of y coords equals #rows
nx=size(X,2); % # of x coords equals #columns
U100=zeros(size(X));
gz100=zeros(size(X));

for j=1:nx
    for i=1:ny
        U100_1=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_1,m_1,G);
        U100_2=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_2,m_2,G);
        U100_3=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_3,m_3,G);
        U100_4=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_4,m_4,G);
        U100_5=grav_pot_point(([X(i,j),Y(i,j),Z100]),xm_5,m_5,G);
        U100(i,j)=U100_1+U100_2+U100_3+U100_4+U100_5;
        gz100_1=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_1,m_1,G);
        gz100_2=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_2,m_2,G);
        gz100_3=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_3,m_3,G);
        gz100_4=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_4,m_4,G);
        gz100_5=grav_eff_point(([X(i,j),Y(i,j),Z100]),xm_5,m_5,G);
        gz100(i,j)=gz100_1+gz100_2+gz100_3+gz100_4+gz100_5;
    end
end

%% get min and max for U and gz for colorbar limits
Umin=min(min(U100));
Umin=min(Umin,min(min(U0)));

Umax=max(max(U100));
Umax=max(Umax,max(max(U0)));

gzmin=min(min(gz100));
gzmin=min(gzmin,min(min(gz0)));

gzmax=max(max(gz100));
gzmax=max(gzmax,max(max(gz0)));
%% Creating figures contour plots for z=0,10 and 100 for meshgrid spacing 25

figure
subplot(3,2,1);
contourf(X,Y,U0);
hold on; % overlay grid points

plot(X,Y, 'xk', 'MarkerSize',2); %creating plot
axis equal; % plot x and y on equal scales
axis([xmin,xmax,ymin,ymax]); %define axis limits

xlabel('x [m]', 'fontweight','bold'); %x, y labels and title
ylabel('y [m]', 'fontweight','bold'); 
title('Gravitational Potential - depth 0m - 25m spacing - 5 random masses (set3)')

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
title('Gravitational Effect - depth 0m - 25m spacing - 5 random masses (set3)')

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
title('Gravitational Potential - depth 10m - 25m spacing - 5 random masses (set3)')

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
title('Gravitational Effect - depth 10m - 25m spacing - 5 random masses (set3)')
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
title('Gravitational Potential - depth 100m - 25m spacing - 5 random masses (set3)')
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
title('Gravitational E - depth 0m - 25m spacing - 5 random masses (set3)')
h_c=colorbar; 
ylabel(h_c, '[m/s^2]','fontweight','bold');
caxis([gzmin,gzmax])

cmap=colormap('parula');
cmap=flipud(cmap);
colormap(cmap);
