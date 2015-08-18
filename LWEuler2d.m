%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lax-Wendroff method to solve 2-D Euler equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Following the ideas of:
% 1. E.F. Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics 
% Manchester U.K., Springer Editorial, 2nd Ed. Chapert 11.
%
% coded by Manuel Diaz, 2012.12.25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; %close all; clc;

%% Parameters
CFL   = 0.20;	% CFL number
tEnd  = 0.10;	% Final time
nE    = 100;    % Number of cells/Elements
n     = 5;      % Number of degrees of freedom
alpha = 0.55;	% Parameter for artificial viscosity
IC    = 'test';	% 19 IC cases are available

% Ratio of specific heats for ideal di-atomic gas
gamma=(n+2)/n;

% Discretize spatial domain
a=0; b=1; dx=(b-a)/nE; nx=nE+1;
c=0; d=1; dy=(d-c)/nE; ny=nE+1;
[x,y] = meshgrid(linspace(a,b,nx),linspace(c,d,ny));

% Set IC
[r0,u0,v0,p0] = Euler_IC2d(x,y,IC);
E0 = p0./((gamma-1)*r0)+0.5*(u0.^2+v0.^2);  % Total Energy
a0 = sqrt(gamma*p0./r0);	% Speed of sound

% Discretize time domain
lambda1 = abs(u0)+a0;
lambda2 = abs(v0)+a0;
dt=CFL*dx/max([lambda1(:);lambda2(:)]);  % using the system's largest eigenvalue
t = 0:dt:tEnd; 

%% Solver Loop
% Load initial condition
r=r0; u=u0; v=v0; p=p0; E=E0; it=0;

% Allocate auxiliary variables
zero=zeros(ny-1,nx-1); one=ones(ny-1,nx-1); 
viscX=zeros(ny-1,nx-1,4); viscY=zeros(ny-1,nx-1,4); qs=zeros(ny-1,nx-1,4);

for tsteps=t
    % Iteration counter
    it=it+1;
    
    % Define vectors q, F and G for every {x(i),y(i)}
    q = cat(3, r, r.*u, r.*v, r.*E);
    F = cat(3, r.*u, r.*u.^2+p, r.*u.*v, u.*(r.*E+p));
    G = cat(3, r.*v, r.*u.*v, r.*v.^2+p, v.*(r.*E+p));
    
    % Compute q* matrix and flow parameters
    qs = 0.25*(q(1:ny-1,2:nx,:) + q(1:ny-1,1:nx-1,:)) ...
        +0.25*(q(2:ny,1:nx-1,:) + q(1:ny-1,1:nx-1,:)) ...
        -dt/(2*dx)*(F(1:ny-1,2:nx,:) - F(1:ny-1,1:nx-1,:)) ...
        -dt/(2*dy)*(G(2:ny,1:nx-1,:) - G(1:ny-1,1:nx-1,:));
    
	% Compute flow properties
    rs=qs(:,:,1);       us=qs(:,:,2)./rs; 
    vs=qs(:,:,3)./rs;   Es=qs(:,:,4)./rs;
    ps=(gamma-1)*rs.*(Es-0.5*(us.^2+vs.^2));
    
    Fs = cat(3, rs.*us, rs.*us.^2+ps, rs.*us.*vs, us.*(rs.*Es+ps));
    Gs = cat(3, rs.*vs, rs.*us.*vs, rs.*vs.^2+ps, vs.*(rs.*Es+ps));
    mx = cat(3,zero,one,one,us);
    
    % Calculate arfiticial viscosity
    for j=1:4
        viscX(:,:,j) = zero; %alpha*dx^2*r(1:nx-1)...
                        %.*abs((u(:,2:nx)-u(:,1:nx-1))/dx)...
                        %.*((u(2:nx)-u(1:nx-1))/dx).*mx(:,:,j);
        viscY(:,:,j) = zero; %alpha*dy^2*r(1:nx-1)...
                        %.*abs((u(2:nx,:)-u(1:nx-1,:))/dy)...
                        %.*((u(2:nx,:)-u(1:nx-1,:))/dy).*mx(:,:,j);
    end
	Fs=Fs-viscX; Gs=Gs-viscY;
    
    % Update F* and q matrix
    q(2:ny-1,2:nx-1,:) = q(2:ny-1,2:nx-1,:) ...
        -dt/dx*(Fs(1:ny-2,2:nx-1,:) - Fs(1:ny-2,1:nx-2,:)) ...
        -dt/dy*(Gs(2:ny-1,1:nx-2,:) - Gs(1:ny-2,1:nx-2,:));
    
    % Neumann BCs
    q(:,1,:)  = q(:,2,:);       q(1,:,:)  = q(2,:,:);
    q(:,nx,:) = q(:,nx-1,:);    q(ny,:,:) = q(ny-1,:,:);
                
	% Compute flow properties
    r=q(:,:,1);   u=q(:,:,2)./r;    v=q(:,:,3)./r;    E=q(:,:,4)./r;
    p=(gamma-1)*r.*(E-0.5*(u.^2+v.^2));
end

%% Calculation of flow parameters
a = sqrt(gamma*p./r);
Mx = u./a;
My = v./a;
U = sqrt(u.^2+v.^2);
M = U./a;
p_ref = 101325;     % Reference air pressure (N/m^2)
rho_ref= 1.225;     % Reference air density (kg/m^3)
s = 1/(gamma-1)*(log(p/p_ref)+gamma*log(rho_ref./r)); 
                    % Entropy w.r.t reference condition
ss = log(p./r.^gamma);
                    % Dimensionless Entropy
r_x = r.*u;         % Mass Flow rate per unit area
r_y = r.*v;         % Mass Flow rate per unit area
e = p./((gamma-1)*r); % internal Energy

%% Final plot
offset=0.05; n=22; % contour lines
s1=subplot(2,3,1); contour(x,y,r,n); axis('square'); xlabel('x(m)'); ylabel('Density (kg/m^3)');
s2=subplot(2,3,2); contour(x,y,U,n); axis('square'); xlabel('x(m)'); ylabel('Velocity Magnitud (m/s)');
s3=subplot(2,3,3); contour(x,y,p,n); axis('square'); xlabel('x(m)'); ylabel('Pressure (Pa)');
s4=subplot(2,3,4); contour(x,y,ss,n);axis('square'); xlabel('x(m)'); ylabel('Entropy/R gas');
s5=subplot(2,3,5); contour(x,y,M,n); axis('square'); xlabel('x(m)'); ylabel('Mach number');
s6=subplot(2,3,6); contour(x,y,e,n); axis('square'); xlabel('x(m)'); ylabel('Internal Energy (kg/m^2s)');
title(s1,'Lax-Wendroff 2-D Euler Solver');