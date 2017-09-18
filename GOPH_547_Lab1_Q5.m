%GOPH 547 - Lab 1 - Forward Modelling for Gravity Anomalies 
%Safian Omar Qureshi
%ID: 10086638
%TA: Ye Sun
%Collaberated with Tian Yu, Ryan Younghwan Ok

%% Question 5a:
clear; clc;

load('anomaly_data.mat');

% finding j parameters
jz_8=find(z(1,1,:)==-8,1); %when z=-8 (z axis positive downward) 
jz_10=find(z(1,1,:)==-10,1); %when z=-10
jz_12=find(z(1,1,:)==-12,1); %when z=-12

jy_neg8=find(y(:,1,1)==-8,1); %when y=-8
jy_0=find(y(:,1,1)==0,1);%when y=0
jy_8=find(y(:,1,1)==8,1);%when y=8

jx_neg16=find(x(1,:,1)==-16,1);%when x=-16
jx_0=find(x(1,:,1)==0,1);%when x=0
jx_16=find(x(1,:,1)==16,1);%when x=16

%load x,y and z grid to plot corresponding rho values in cross-section.
ny=size(x,1); % number of y points is the number of ROWS
nx=size(x,2); % number of x points is the number of COLUMNS
nz=size(z,3); % number of z points is the number of 'PAGES' (3rd 'axis' dimension)

% For z=-8
x_cur_z8=zeros(nx,ny); x_cur_z8(:,:)=x(:,:,jz_8);
y_cur_z8=zeros(nx,ny); y_cur_z8(:,:)=y(:,:,jz_8);
rho_cur_z8=zeros(nx,ny); rho_cur_z8(:,:)=rho(:,:,jz_8);

% For z=-10
x_cur_z10=zeros(nx,ny); x_cur_z10(:,:)=x(:,:,jz_10);
y_cur_z10=zeros(nx,ny); y_cur_z10(:,:)=y(:,:,jz_10);
rho_cur_z10=zeros(nx,ny); rho_cur_z10(:,:)=rho(:,:,jz_10);

% For z=-12
x_cur_z12=zeros(nx,ny); x_cur_z12(:,:)=x(:,:,jz_12);
y_cur_z12=zeros(nx,ny); y_cur_z12(:,:)=y(:,:,jz_12);
rho_cur_z12=zeros(nx,ny); rho_cur_z12(:,:)=rho(:,:,jz_12);

% For y=-8
x_cur_yneg8=zeros(nx,nz); x_cur_yneg8(:,:)=x(jy_neg8,:,:);
z_cur_yneg8=zeros(nx,nz); z_cur_yneg8(:,:)=z(jy_neg8,:,:);
rho_cur_yneg8=zeros(nx,nz); rho_cur_yneg8(:,:)=rho(jy_neg8,:,:);

% For y=0
x_cur_y0=zeros(nx,nz); x_cur_y0(:,:)=x(jy_0,:,:);
z_cur_y0=zeros(nx,nz); z_cur_y0(:,:)=z(jy_0,:,:);
rho_cur_y0=zeros(nx,nz); rho_cur_y0(:,:)=rho(jy_0,:,:);

% For y=8
x_cur_y8=zeros(nx,nz); x_cur_y8(:,:)=x(jy_8,:,:);
z_cur_y8=zeros(nx,nz); z_cur_y8(:,:)=z(jy_8,:,:);
rho_cur_y8=zeros(nx,nz); rho_cur_y8(:,:)=rho(jy_8,:,:);

% For x=-16
y_cur_xneg16=zeros(ny,nz); y_cur_xneg16(:,:)=y(:,jx_neg16,:);
z_cur_xneg16=zeros(ny,nz); z_cur_xneg16(:,:)=z(:,jx_neg16,:);
rho_cur_xneg16=zeros(ny,nz); rho_cur_xneg16(:,:)=rho(:,jx_neg16,:);

% For x=0
y_cur_x0=zeros(ny,nz); y_cur_x0(:,:)=y(:,jx_0,:);
z_cur_x0=zeros(ny,nz); z_cur_x0(:,:)=z(:,jx_0,:);
rho_cur_x0=zeros(ny,nz); rho_cur_x0(:,:)=rho(:,jx_0,:);

% For x=16
y_cur_x16=zeros(ny,nz); y_cur_x16(:,:)=y(:,jx_16,:);
z_cur_x16=zeros(ny,nz); z_cur_x16(:,:)=z(:,jx_16,:);
rho_cur_x16=zeros(ny,nz); rho_cur_x16(:,:)=rho(:,jx_16,:);

% plots the assigned sections
figure; % For z=-8
subplot(3,3,1)
contourf(x_cur_z8,y_cur_z8,rho_cur_z8);
h_c = colorbar; % add a colorbar
cmap=colormap('parula');
cmap=flipud(cmap);
colormap(cmap);
axis equal;
axis([-20,20,-12,12]);
xlabel('x [m]', 'fontweight','bold'); % add labels
ylabel('y [m]', 'fontweight','bold');
title('Density anomaly map z = -8m', 'fontweight','bold');
ylabel(h_c, '\rho [g/cm^3]', 'fontweight','bold');
%%
subplot(3,3,2) % For z=-10
contourf(x_cur_z10,y_cur_z10,rho_cur_z10);
h_c = colorbar; % add a colorbar
cmap=colormap('parula');
cmap=flipud(cmap);
colormap(cmap);
axis equal;
axis([-20,20,-12,12]);
xlabel('x [m]', 'fontweight','bold'); % add labels
ylabel('y [m]', 'fontweight','bold');
title('Density anomaly map z = -10m', 'fontweight','bold');
ylabel(h_c, '\rho [g/cm^3]', 'fontweight','bold');

subplot(3,3,3) % For z=-12
contourf(x_cur_z12,y_cur_z12,rho_cur_z12);
h_c = colorbar; % add a colorbar
cmap=colormap('parula');
cmap=flipud(cmap);
colormap(cmap);
axis equal;
axis([-20,20,-12,12]);
xlabel('x [m]', 'fontweight','bold'); % add labels
ylabel('y [m]', 'fontweight','bold');
title('Density anomaly map z = -12 m', 'fontweight','bold');
ylabel(h_c, '\rho [g/cm^3]', 'fontweight','bold');

subplot(3,3,4)% For y=-8
contourf(x_cur_yneg8,z_cur_yneg8,rho_cur_yneg8);
h_c = colorbar; % add a colorbar
cmap=colormap('parula');
cmap=flipud(cmap);
colormap(cmap);
axis equal;
axis([-20,20,-15,-5]);
xlabel('x [m]', 'fontweight','bold'); % add labels
ylabel('z [m]', 'fontweight','bold');
title('Density anomaly map y = -8m', 'fontweight','bold');
ylabel(h_c, '\rho [g/cm^3]', 'fontweight','bold');

subplot(3,3,5)% For y=0
contourf(x_cur_y0,z_cur_y0,rho_cur_y0);
h_c = colorbar; % add a colorbar
cmap=colormap('parula');
cmap=flipud(cmap);
colormap(cmap);
axis equal;
axis([-20,20,-15,-5]);
xlabel('x [m]', 'fontweight','bold'); % add labels
ylabel('z [m]', 'fontweight','bold');
title('Density anomaly map y = 0m', 'fontweight','bold');
ylabel(h_c, '\rho [g/cm^3]', 'fontweight','bold');

subplot(3,3,6)% For y=8
contourf(x_cur_y8,z_cur_y8,rho_cur_y8);
h_c = colorbar; % add a colorbar
cmap=colormap('parula');
cmap=flipud(cmap);
colormap(cmap);
axis equal;
axis([-20,20,-15,-5]);
xlabel('x [m]', 'fontweight','bold'); % add labels
ylabel('z [m]', 'fontweight','bold');
title('Density anomaly map y = 8m', 'fontweight','bold');
ylabel(h_c, '\rho [g/cm^3]', 'fontweight','bold');

subplot(3,3,6)% For y=8
contourf(x_cur_y8,z_cur_y8,rho_cur_y8);
h_c = colorbar; % add a colorbar
cmap=colormap('parula');
cmap=flipud(cmap);
colormap(cmap);
axis equal;
axis([-20,20,-15,-5]);
xlabel('x [m]', 'fontweight','bold'); % add labels
ylabel('z [m]', 'fontweight','bold');
title('Density anomaly map y = 8m', 'fontweight','bold');
ylabel(h_c, '\rho [g/cm^3]', 'fontweight','bold');

subplot(3,3,7)% For x=-16
contourf(y_cur_xneg16,z_cur_xneg16,rho_cur_xneg16);
h_c = colorbar; % add a colorbar
cmap=colormap('parula');
cmap=flipud(cmap);
colormap(cmap);
axis equal;
axis([-20,20,-15,-5]);
xlabel('y [m]', 'fontweight','bold'); % add labels
ylabel('z [m]', 'fontweight','bold');
title('Density anomaly map x = -16m', 'fontweight','bold');
ylabel(h_c, '\rho [g/cm^3]', 'fontweight','bold');

subplot(3,3,8)% For x=0
contourf(y_cur_x0,z_cur_x0,rho_cur_x0);
h_c = colorbar; % add a colorbar
cmap=colormap('parula');
cmap=flipud(cmap);
colormap(cmap);
axis equal;
axis([-20,20,-15,-5]);
xlabel('y [m]', 'fontweight','bold'); % add labels
ylabel('z [m]', 'fontweight','bold');
title('Density anomaly map x = 0m', 'fontweight','bold');
ylabel(h_c, '\rho [g/cm^3]', 'fontweight','bold');
%%
subplot(3,3,9)% For x=16
contourf(y_cur_x16,z_cur_x16,rho_cur_x16);
h_c = colorbar; % add a colorbar
cmap=colormap('parula');
cmap=flipud(cmap);
colormap(cmap);
axis equal;
axis([-20,20,-15,-5]);
xlabel('y [m]', 'fontweight','bold'); % add labels
ylabel('z [m]', 'fontweight','bold');
title('Density anomaly map x = 16m', 'fontweight','bold');
ylabel(h_c, '\rho [g/cm^3]', 'fontweight','bold');

%% Question 5b - total mass

for i=1:length(x(:,1,1)) %going over all 3 indices
    for j=1:length(y(1,:,1))
        for k=1:length(z(1,1,:))
            m(i,j,k)=rho(i,j,k)*8*10^3; %using density mass volume relation, applying unit conversion
        end
    end
end

TotM=sum(sum(sum(m)));% summing all indices 
%% Question 5c 
%survey 5 grid spacing

G=6.674*10^(-11);
dx=5; xmin=-100; xmax=100;
dy=5; ymin=-100; ymax=100;

x0=xmin:dx:xmax;
y0=ymin:dy:ymax;

[X,Y]=meshgrid(x0,y0);


Z0=0;
ny=size(X,1); 
nx=size(X,2); 
gz0=zeros(size(X));
ind_nz=find(rho);% finding non-zero density
nnz=length(ind_nz); % # of non-zero position

for j=1:nx
    for i=1:ny
        for k=1:nnz
            kk=ind_nz(k); % index next non-zero
            xm = [ x(kk), y(kk), z(kk) ];% Coordinates of non-zero
            dm = rho(kk)*8*10^3; % v=(2m*2m*2m)=8 m^3, x10^3 to convert g to kg; 
            gz0(i,j)=gz0(i,j)+grav_eff_point(([X(i,j),Y(i,j),Z0]),xm,dm,G);
        end
    end
end



Z100=100;
ny=size(X,1);
nx=size(X,2); 
gz100=zeros(size(X));
ind_nz=find(rho);% finding non-zero density
nnz=length(ind_nz); % # of non-zero loc

for j=1:nx
    for i=1:ny
        for k=1:nnz
            kk=ind_nz(k); % index next non-zero
            xm = [ x(kk), y(kk), z(kk) ];% Coordinates of non-zero
            dm = rho(kk)*8*10^3; %v=(2m*2m*2m)=8 m^3, x10^3 to convert g to kg; 
            gz100(i,j)=gz100(i,j)+grav_eff_point(([X(i,j),Y(i,j),Z100]),xm,dm,G);
        end
    end
end

% min and max gz for colorbar limits
gzmin=min(min(gz100)); 
gzmin=min(gzmin,min(min(gz0)));

gzmax=max(max(gz100));
gzmax=max(gzmax,max(max(gz0)));

% Create figure for z=0,100 and dx/dy=5,25
figure
subplot(2,2,1);
contourf(X,Y,gz0);
hold on; 
plot(X,Y,'xk','MarkerSize',2);
axis equal; 
axis([xmin,xmax,ymin,ymax]); %define axis limits
xlabel('x (m)','fontweight','bold'); 
ylabel('y (m)','fontweight','bold'); 
title('Gravity Effect for a depth of 0m - 5m spacing')
h_c=colorbar; %Add colorbar
ylabel(h_c, '[m/s^2]','fontweight','bold');
caxis([gzmin,gzmax])

subplot(2,2,2);
contourf(X,Y,gz100);
hold on; 
plot(X,Y,'xk','MarkerSize',2);
axis equal; 
axis([xmin,xmax,ymin,ymax]); %define axis limits
xlabel('x (m)','fontweight','bold');
ylabel('y (m)','fontweight','bold'); 
title('Gravity Effect for a depth of 100m - 5m spacing')
h_c=colorbar; 
ylabel(h_c, '[m/s^2]','fontweight','bold');
caxis([gzmin,gzmax])

%Set up meshgrid 25m spacing

dx=25; xmin=-100; xmax=100;
dy=25; ymin=-100; ymax=100;

x0=xmin:dx:xmax;
y0=ymin:dy:ymax;

[X,Y]=meshgrid(x0,y0);

%Loop over all survey points for  z=0 
Z0=0;
ny=size(X,1); 
nx=size(X,2); 
gz0=zeros(size(X));
ind_nz=find(rho);% finding the non-zero density
nnz=length(ind_nz); % # of non-zero position
for j=1:nx
    for i=1:ny
        for k=1:nnz
            kk=ind_nz(k); % index next non-zero
            xm = [ x(kk), y(kk), z(kk) ];% coordinates of non-zero
            dm = rho(kk)*8*10^3; %v=(2m*2m*2m)=8 m^3, x10^3 to convert g to kg; 
            gz0(i,j)=gz0(i,j)+grav_eff_point(([X(i,j),Y(i,j),Z0]),xm,dm,G);
        end
    end
end


%100m case 
Z100=100;
ny=size(X,1);
nx=size(X,2); 
gz100=zeros(size(X));
ind_nz=find(rho);
nnz=length(ind_nz); 

for j=1:nx
    for i=1:ny
        for k=1:nnz
            kk=ind_nz(k); % index next non-zero
            xm = [ x(kk), y(kk), z(kk) ];% Coordinates of non-zero
            dm = rho(kk)*8*10^3; %v=(2m*2m*2m)=8 m^3, x10^3 to convert g to kg; 
            gz100(i,j)=gz100(i,j)+grav_eff_point(([X(i,j),Y(i,j),Z100]),xm,dm,G);
        end
    end
end


% get min and max gz for colorbar limits
gzmin=min(min(gz100));
gzmin=min(gzmin,min(min(gz0)));

gzmax=max(max(gz100));
gzmax=max(gzmax,max(max(gz0)));



subplot(2,2,3);
contourf(X,Y,gz0);
hold on; %overlay grid pts
plot(X,Y,'xk','MarkerSize',2);
axis equal; %plot x & y on equal scales
axis([xmin,xmax,ymin,ymax]); %define axis limits
xlabel('x (m)','fontweight','bold'); %Adding axis labels/title
ylabel('y (m)','fontweight','bold'); 
title('Gravity Effect for a depth of 0m - 25m spacing')
h_c=colorbar; %Adding colorbar
ylabel(h_c, '[m/s^2]','fontweight','bold');
caxis([gzmin,gzmax])

subplot(2,2,4);
contourf(X,Y,gz100);
hold on; %overlay grid pts
plot(X,Y,'xk','MarkerSize',2);
axis equal; %plot x & y on equal scales
axis([xmin,xmax,ymin,ymax]); %define axis limits
xlabel('x (m)','fontweight','bold'); %Adding axis labels/title
ylabel('y (m)','fontweight','bold'); 
title('Gravity Effect for a depth of 100m - 25m spacing')
h_c=colorbar; %Adding colorbar
ylabel(h_c, '[m/s^2]','fontweight','bold');
caxis([gzmin,gzmax])
cmap=colormap('parula');
cmap=flipud(cmap);
colormap(cmap);

%% Part D
%Setting up survey meshgrid for 5m spacing
G=6.674*10^(-11);
dx=5; xmin=-100; xmax=100;
dy=5; ymin=-100; ymax=100;

x0=xmin:dx:xmax;
y0=ymin:dy:ymax;

[X,Y]=meshgrid(x0,y0);

%Loop over survey pts for the z=0 case. 
Z0=0;
ny=size(X,1); % # of y coords = # of rows
nx=size(X,2); % # of x coords = # of columns
gz0=zeros(size(X));
ind_nz=find(rho);% finding the non-zero density
nnz=length(ind_nz); % # of non-zero loc.
for j=1:nx
    for i=1:ny
        for k=1:nnz
            kk=ind_nz(k); % indext next non-zero
            xm = [ x(kk), y(kk), z(kk) ];% Coordinates of non-zero
            dm = rho(kk)*8*10^3; %v=(2m*2m*2m)=8 m^3, x10^3 to convert g to kg; 
            gz0(i,j)=gz0(i,j)+grav_eff_point(([X(i,j),Y(i,j),Z0]),xm,dm,G);
        end
    end
end


%Loop over survey pts for the z=100 case. 
Z100=100;
ny=size(X,1); % # of y coords = # of rows
nx=size(X,2); % # of x coords = # of columns
gz100=zeros(size(X));
ind_nz=find(rho);% finding the non-zero density
nnz=length(ind_nz); % # of non-zero loc.
for j=1:nx
    for i=1:ny
        for k=1:nnz
            kk=ind_nz(k); % indext next non-zero
            xm = [ x(kk), y(kk), z(kk) ];% Coordinates of non-zero
            dm = rho(kk)*8*10^3; %v=(2m*2m*2m)=8 m^3, x10^3 to convert g to kg; 
            gz100(i,j)=gz100(i,j)+grav_eff_point(([X(i,j),Y(i,j),Z100]),xm,dm,G);
        end
    end
end

%Loop over survey pts for the z=1 case. 
z1=1;
ny=size(X,1); % # of y coords = # of rows
nx=size(X,2); % # of x coords = # of columns
gz1=zeros(size(X));
ind_nz=find(rho);% finding the non-zero density
nnz=length(ind_nz); % # of non-zero loc.
for j=1:nx
    for i=1:ny
        for k=1:nnz
            kk=ind_nz(k); % indext next non-zero
            xm = [ x(kk), y(kk), z(kk) ];% Coordinates of non-zero
            dm = rho(kk)*8*10^3; %v=(2m*2m*2m)=8 m^3, x10^3 to convert g to kg; 
            gz1(i,j)=gz1(i,j)+grav_eff_point(([X(i,j),Y(i,j),z1]),xm,dm,G);
        end
    end
end


%Loop over survey pts for the z=110 case. 
z110=110;
ny=size(X,1); % # of y coords = # of rows
nx=size(X,2); % # of x coords = # of columns
gz110=zeros(size(X));
ind_nz=find(rho);% finding the non-zero density
nnz=length(ind_nz); % # of non-zero loc.
for j=1:nx
    for i=1:ny
        for k=1:nnz
            kk=ind_nz(k); % indext next non-zero
            xm = [ x(kk), y(kk), z(kk) ];% Coordinates of non-zero
            dm = rho(kk)*8*10^3; %v=(2m*2m*2m)=8 m^3, x10^3 to convert g to kg; 
            gz110(i,j)=gz110(i,j)+grav_eff_point(([X(i,j),Y(i,j),z110]),xm,dm,G);
        end
    end
end


% Derivative

dz1=1;
dg1=(-gz1+gz0);
dz2=10;
dg2=(-gz110+gz100)/dz2;


% Get min/max gz for colorbar limits
dgmin=min(min(dg1));
dgmin=min(dgmin,min(min(dg2)));
dgmax=max(max(dg1));
dgmax=max(dgmax,max(max(dg2)));

figure
subplot(2,2,1);
contourf(X,Y,dg1);
hold on; %overlay grid pts
plot(X,Y,'xk','MarkerSize',2);
axis equal; %plot x & y on equal scales
axis([xmin,xmax,ymin,ymax]); %define axis limits
xlabel('x (m)','fontweight','bold'); %Adding axis labels/title
ylabel('y (m)','fontweight','bold'); 
title('First Partial Derivative between 0 and 1m - 5m spacing')
h_c=colorbar; %Adding colorbar
ylabel(h_c, '[m/s^2]','fontweight','bold');
caxis([dgmin,dgmax])

subplot(2,2,2);
contourf(X,Y,dg2);
hold on; %overlay grid pts
plot(X,Y,'xk','MarkerSize',2);
axis equal; %plot x & y on equal scales
axis([xmin,xmax,ymin,ymax]); %define axis limits
xlabel('x (m)','fontweight','bold'); %Adding axis labels/title
ylabel('y (m)','fontweight','bold'); 
title('First Partial Derivative between 100 and 110m - 5m spacing')
h_c=colorbar; %Adding colorbar
ylabel(h_c, '[m/s^2]','fontweight','bold');
caxis([dgmin,dgmax])

%Setting up survey meshgrid for 25m spacing
G=6.674*10^(-11);
dx=25; xmin=-100; xmax=100;
dy=25; ymin=-100; ymax=100;

x0=xmin:dx:xmax;
y0=ymin:dy:ymax;

[X,Y]=meshgrid(x0,y0);

%loop over survey pts for the z=0 case. 
Z0=0;
ny=size(X,1);
nx=size(X,2); 
gz0=zeros(size(X));
ind_nz=find(rho);% finding  non-zero density
nnz=length(ind_nz); % # of non-zero positions
for j=1:nx
    for i=1:ny
        for k=1:nnz
            kk=ind_nz(k); % indext next non-zero
            xm = [ x(kk), y(kk), z(kk) ];% Coordinates of non-zero
            dm = rho(kk)*8*10^3; %v=(2m*2m*2m)=8 m^3, x10^3 to convert g to kg; 
            gz0(i,j)=gz0(i,j)+grav_eff_point(([X(i,j),Y(i,j),Z0]),xm,dm,G);
        end
    end
end


%Loop over survey pts for the z=100 case. 

Z100=100;
ny=size(X,1); % # of y coords = # of rows
nx=size(X,2); % # of x coords = # of columns
gz100=zeros(size(X));
ind_nz=find(rho);% finding the non-zero density
nnz=length(ind_nz); % # of non-zero loc.
for j=1:nx
    for i=1:ny
        for k=1:nnz
            kk=ind_nz(k); % indext next non-zero
            xm = [ x(kk), y(kk), z(kk) ];% Coordinates of non-zero
            dm = rho(kk)*8*10^3; %v=(2m*2m*2m)=8 m^3, x10^3 to convert g to kg; 
            gz100(i,j)=gz100(i,j)+grav_eff_point(([X(i,j),Y(i,j),Z100]),xm,dm,G);
        end
    end
end

%Loop over survey pts for the z=1 case. 
z1=1;
ny=size(X,1); % # of y coords = # of rows
nx=size(X,2); % # of x coords = # of columns
gz1=zeros(size(X));
ind_nz=find(rho);% finding the non-zero density
nnz=length(ind_nz); % # of non-zero loc.
for j=1:nx
    for i=1:ny
        for k=1:nnz
            kk=ind_nz(k); % indext next non-zero
            xm = [ x(kk), y(kk), z(kk) ];% Coordinates of non-zero
            dm = rho(kk)*8*10^3; %v=(2m*2m*2m)=8 m^3, x10^3 to convert g to kg; 
            gz1(i,j)=gz1(i,j)+grav_eff_point(([X(i,j),Y(i,j),z1]),xm,dm,G);
        end
    end
end


%Loop over survey pts for the z=110 case. 
z110=110;
ny=size(X,1); % # of y coords = # of rows
nx=size(X,2); % # of x coords = # of columns
gz110=zeros(size(X));
ind_nz=find(rho);% finding the non-zero density
nnz=length(ind_nz); % # of non-zero loc.

for j=1:nx
    for i=1:ny
        for k=1:nnz
            kk=ind_nz(k); % indext next non-zero
            xm = [ x(kk), y(kk), z(kk) ];% Coordinates of non-zero
            dm = rho(kk)*8*10^3; %v=(2m*2m*2m)=8 m^3, x10^3 to convert g to kg; 
            gz110(i,j)=gz110(i,j)+grav_eff_point(([X(i,j),Y(i,j),z110]),xm,dm,G);
        end
    end
end


% Derivative

dz1=1;
dg1=(-gz1+gz0);
dz2=10;
dg2=(-gz110+gz100)/dz2;


% Get min/max gz for colorbar limits
dgmin=min(min(dg1));
dgmin=min(dgmin,min(min(dg2)));
dgmax=max(max(dg1));
dgmax=max(dgmax,max(max(dg2)));

subplot(2,2,3);
contourf(X,Y,dg1);
hold on; %overlay grid pts
plot(X,Y,'xk','MarkerSize',2);
axis equal; %plot x & y on equal scales
axis([xmin,xmax,ymin,ymax]); %define axis limits
xlabel('x (m)','fontweight','bold'); %Adding axis labels/title
ylabel('y (m)','fontweight','bold'); 
title('First Partial Derivative between 0 and 1m - 25m spacing')
h_c=colorbar; %Adding colorbar
ylabel(h_c, '[m/s^2]','fontweight','bold');
caxis([dgmin,dgmax])

subplot(2,2,4);
contourf(X,Y,dg2);
hold on; %overlay grid pts
plot(X,Y,'xk','MarkerSize',2);
axis equal; %plot x & y on equal scales
axis([xmin,xmax,ymin,ymax]); %define axis limits
xlabel('x (m)','fontweight','bold'); %Adding axis labels/title
ylabel('y (m)','fontweight','bold'); 
title('First Partial Derivative between 100 and 110m - 25m spacing')
h_c=colorbar; %Adding colorbar
ylabel(h_c, '[m/s^2]','fontweight','bold');
caxis([dgmin,dgmax])
cmap=colormap('parula');
cmap=flipud(cmap);
colormap(cmap);
