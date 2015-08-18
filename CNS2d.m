%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 A 2-D Compresssible Navier-Stokes Equations
%                      by Manuel Diaz, NTU, 20.03.2015
%
% A finite volume implementaion for the compresssible navier-stokes (CNS)
% equations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% [1] Are Skoinen, Cartesian Grid methods for the compressible
%     Navier-Stokes equations. Master Thesis, NTNU Norway, Spring 2012.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% Method Parameters
CFL     = 0.60;     % CFL number
VNN     = 0.40;     % Von Neumann number
tEnd    = 0.15;     % Final time
nx      = 100;      % Number of cells/Elements in x
ny      = 002;      % Number of cells/Elements in y
n       = 5;        % Number of degrees of freedom
IC      ='Sod';     % 'Sod', 'sShock'
limiter ='MC';      % MM, MC, VA.
fluxMth ='RUS';     % LF, RUS, ROE, HLLE, RHLLE(!).
plot_fig= 1;        % 1:visualize evolution 

%% Physical parameters
 gamma = (n+2)/n; % Ratio of specific heats for ideal di-atomic gas
mu_ref = 1.716E-5;  %[]
 T_ref = 273.1;     %[K]
     R = 287;       %[N.m/kg]
    Pr = 0.72;      %[-] Prandtl number
    cp = R*gamma/(gamma-1);
    %mu = mu_ref*(T_ref+110)/(T+110)*(T/T_ref)^1.5; % viscosity \mu
    % k = mu*cp/Pr; % thermal conductivity
    mu = 0.01;
     k = 0.01;

%% Preprocess
% Build a (NCFV) in a rectangular domain of size [-Lx,Lx]x[-Ly,Ly] using a
% cartesian grid of [NxM] points with to distribute [(nx)x(ny)] cells of
% equidistant spacing dx = 2*Lx/(nx) and dy = 2*Ly/(ny)
Lx=1; dx=2*Lx/nx; xc=-Lx+dx/2:dx:Lx;
Ly=1; dy=2*Ly/ny; yc=-Ly+dy/2:dy:Ly;
[x,y] = meshgrid(xc,yc); % cell centered values

%% Load IC
[r0,u0,v0,p0] = CNS_IC2d(x,y,gamma,IC);

% Temperature, Total Energy density, Enthalpy density:
U0 = sqrt(u0.^2+v0.^2);
E0 = p0./((gamma-1)*r0)+0.5*(u0.^2+v0.^2);  % Total Energy
c0 = sqrt(gamma*p0./r0);                    % Speed of sound
Q0 = cat(3, r0, r0.*u0, r0.*v0, r0.*E0);    % initial state

% Set q-array & adjust grid for ghost cells
nx=nx+2; ny=ny+2; q0=zeros(ny,nx,4); q0(2:ny-1,2:nx-1,1:4)=Q0;

% Boundary Conditions in ghost cells
q0(:,1,:)=q0(:,2,:); q0(:,nx,:)=q0(:,nx-1,:);   % Natural BCs
q0(1,:,:)=q0(2,:,:); q0(ny,:,:)=q0(ny-1,:,:);   % Natural BCs

% Eigenvalues of the system
lambda1=abs(U0)+c0; lambda3=abs(U0)-c0;
smax = max([lambda1(:),lambda3(:)]);

% Discretize time domain
lambda1=U0+c0; lambda2=U0-c0; a0=max(abs([lambda1(:);lambda2(:)])); 
dtInvicid = CFL*min(dx./a0,dy./a0);
lambda3=4/3*mu./r0; lambda4=gamma*mu./(r0*Pr); b0=max([lambda3(:);lambda4(:)]);
dtViscous = VNN/((1/dx^2+1/dy^2)*b0);
dt0=min(dtInvicid,dtViscous);

%% Solver Loop
% load initial conditions
q=q0; it=0; t=0; dt=dt0; a=a0;

% solve loop
tic
while t < tEnd
   % RK3 step 1
   qs = q - dt*MUSCL_EulerSys2d(q,a,gamma,dx,dy,nx,ny,limiter,fluxMth) ...
           + dt*fluxCNS2d(q,gamma,dx,dy,nx,ny,mu,k);
   
   qs(:,1,:)=qs(:,2,:); qs(:,nx,:)=qs(:,nx-1,:);	% Natural BCs
   qs(1,:,:)=qs(2,:,:); qs(ny,:,:)=qs(ny-1,:,:);	% Natural BCs
   
   % RK3 step 2
   qs2 = (3*q + qs - dt*MUSCL_EulerSys2d(qs,a,gamma,dx,dy,nx,ny,limiter,fluxMth) ...
           + dt*fluxCNS2d(qs,gamma,dx,dy,nx,ny,mu,k))/4;
   
   qs2(:,1,:)=qs2(:,2,:); qs2(:,nx,:)=qs2(:,nx-1,:);	% Natural BCs
   qs2(1,:,:)=qs2(2,:,:); qs2(ny,:,:)=qs2(ny-1,:,:);    % Natural BCs
   
   % RK3 step 3 / update solution
   q = (q + 2*(qs2 - dt*MUSCL_EulerSys2d(qs2,a,gamma,dx,dy,nx,ny,limiter,fluxMth) ...
           + dt*fluxCNS2d(qs2,gamma,dx,dy,nx,ny,mu,k)))/3;
   
   q(:,1,:)=q(:,2,:); q(:,nx,:)=q(:,nx-1,:);	% Natural BCs
   q(1,:,:)=q(2,:,:); q(ny,:,:)=q(ny-1,:,:);	% Natural BCs
   
   % compute flow properties
    r=q(:,:,1); u=q(:,:,2)./r; v=q(:,:,3)./r; E=q(:,:,4)./r;
    p=(gamma-1)*r.*(E-0.5*(u.^2+v.^2)); c=sqrt(gamma*p./r);
   
   % Update dt and time
    U=sqrt(u.^2+v.^2); lambda1=U+c; lambda2=U-c;
    a = max(abs([lambda1(:);lambda2(:)]));
    lambda3=4/3*mu./r0; lambda4=gamma*mu./(r0*Pr); 
    b = max([lambda3(:);lambda4(:)]);
    dt=min(CFL*min(dx/a,dy/a),VNN/((1/dx^2+1/dy^2)*b)); 
    if t+dt>tEnd; dt=tEnd-t; end; t=t+dt; it=it+1;
   
   % plot evolution
    if rem(it,1) == 0
        if plot_fig == 1;
            %subplot(2,2,1); contourf(x,y,r(2:ny-1,2:nx-1)); axis('square');
            %subplot(2,2,2); contourf(x,y,u(2:ny-1,2:nx-1)); axis('square');
            %subplot(2,2,3); contourf(x,y,v(2:ny-1,2:nx-1)); axis('square');
            %subplot(2,2,4); contourf(x,y,p(2:ny-1,2:nx-1)); axis('square');
            subplot(2,2,1); plot(x(1,:),r(1,2:nx-1),'.b');
            subplot(2,2,2); plot(x(1,:),u(1,2:nx-1),'.m'); 
            subplot(2,2,3); plot(x(1,:),p(1,2:nx-1),'.k'); 
            subplot(2,2,4); plot(x(1,:),E(1,2:nx-1),'.r');
            drawnow
        end
    end
   
end
cputime = toc;

% Remove ghost cells
q=q(2:ny-1,2:nx-1,1:4); nx=nx-2; ny=ny-2; 

% compute flow properties
r=q(:,:,1); u=q(:,:,2)./r; v=q(:,:,3)./r; E=q(:,:,4)./r; p=(gamma-1)*r.*(E-0.5*(u.^2+v.^2));

%% 
fprintf('completed!\n\n');