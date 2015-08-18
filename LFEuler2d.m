%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lax-Friedrichs method to solve 2-D Euler equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Following the ideas of:
% [1] E.F. Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics 
%     Manchester U.K., Springer Editorial, 2nd Ed., 1999. Chapert 11.
% [2] Randall J. Leveque, Finite Volume Method for Hyperbolic Problems.,
%     Cambridge University Press. 2nd Ed., 2004. Chapter 4.
%
% coded by Manuel Diaz, 2012.12.25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; %close all; clc;

%% Parameters
CFL     = 0.50;	% CFL number
tEnd    = 0.30;	% Final time
nE      = 400;  % Number of cells/Elements
n       = 5;	% Number of degrees of freedom
IC      = 03;   % 19 IC cases are available
plot_fig= 0;    % 1:visualize evolution 

% Ratio of specific heats for ideal di-atomic gas
gamma=(n+2)/n;

% Discretize spatial domain
a=0; b=1; dx=(b-a)/nE; nx=nE+1;
c=0; d=1; dy=(d-c)/nE; ny=nE+1;
[x,y] = meshgrid(linspace(a,b,nx),linspace(c,d,ny));

% Set IC
[r0,u0,v0,p0] = Euler_IC2d(x,y,IC);
E0 = p0./((gamma-1)*r0)+0.5*(u0.^2+v0.^2);  % Total Energy
a0 = sqrt(gamma*p0./r0);            % Speed of sound

% Discretize time domain
lambda1 = abs(u0)+a0;
lambda2 = abs(v0)+a0;
dt=CFL*dx/max([lambda1(:);lambda2(:)]);  % using the system's largest eigenvalue
t = 0:dt:tEnd; 

%% Solver Loop
% Load initial condition
r=r0; u=u0; v=v0; p=p0; E=E0; it=0; 

for tsteps=t
    % iteration counter
    it=it+1;
    
    % define vectors q, F and G for every {x(i),y(i)}
    q = cat(3, r, r.*u, r.*v, r.*E);
    F = cat(3, r.*u, r.*u.^2+p, r.*u.*v, u.*(r.*E+p));
    G = cat(3, r.*v, r.*u.*v, r.*v.^2+p, v.*(r.*E+p));
    
    % update q matrix and flow parameters
    q(2:ny-1,2:nx-1,:) = 0.25*(q(2:ny-1,3:nx,:) + q(2:ny-1,1:nx-2,:)) ...
                       + 0.25*(q(3:ny,2:nx-1,:) + q(1:ny-2,2:nx-1,:)) ...
                    -dt/(2*dx)*(F(2:ny-1,3:nx,:) - F(2:ny-1,1:nx-2,:)) ...
                    -dt/(2*dy)*(G(3:ny,2:nx-1,:) - G(1:ny-2,2:nx-1,:));
    
	% Neumann BCs
    q(:,1,:)  = q(:,2,:);       q(1,:,:)  = q(2,:,:);
    q(:,nx,:) = q(:,nx-1,:);    q(ny,:,:) = q(ny-1,:,:);
                
	% compute flow properties
    r=q(:,:,1); u=q(:,:,2)./r; v=q(:,:,3)./r; E=q(:,:,4)./r;
    p=(gamma-1)*r.*(E-0.5*(u.^2+v.^2));
    
    % Plot figure
    if rem(it,10) == 0
        if plot_fig == 1;
            subplot(2,2,1); contourf(x,y,r); axis('square');
            subplot(2,2,2); contourf(x,y,u); axis('square');
            subplot(2,2,3); contourf(x,y,v); axis('square');
            subplot(2,2,4); contourf(x,y,p); axis('square');
            drawnow
        end
    end
end

%% Calculation of flow parameters
a = sqrt(gamma*p./r); Mx = u./a; My = v./a; U = sqrt(u.^2+v.^2); M = U./a;
p_ref = 101325;         % Reference air pressure (N/m^2)
rho_ref= 1.225;         % Reference air density (kg/m^3)
s = 1/(gamma-1)*(log(p/p_ref)+gamma*log(rho_ref./r)); 
                        % Entropy w.r.t reference condition
ss = log(p./r.^gamma);  % Dimensionless Entropy
r_x = r.*u;             % Mass Flow rate per unit area
r_y = r.*v;             % Mass Flow rate per unit area
e = p./((gamma-1)*r); % internal Energy

%% Final plot
offset=0.05; n=22; % contour lines
s1=subplot(2,3,1); contour(x,y,r,n); axis('square'); xlabel('x(m)'); ylabel('Density (kg/m^3)');
s2=subplot(2,3,2); contour(x,y,U,n); axis('square'); xlabel('x(m)'); ylabel('Velocity Magnitud (m/s)');
s3=subplot(2,3,3); contour(x,y,p,n); axis('square'); xlabel('x(m)'); ylabel('Pressure (Pa)');
s4=subplot(2,3,4); contour(x,y,ss,n);axis('square'); xlabel('x(m)'); ylabel('Entropy/R gas');
s5=subplot(2,3,5); contour(x,y,M,n); axis('square'); xlabel('x(m)'); ylabel('Mach number');
s6=subplot(2,3,6); contour(x,y,e,n); axis('square'); xlabel('x(m)'); ylabel('Internal Energy (kg/m^2s)');
title(s1,'Lax-Friedrichs Euler 2-D Solver');