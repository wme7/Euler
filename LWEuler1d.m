%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lax-Wendroff method to solve 1-D Euler equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Following the ideas of:
% 1. E.F. Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics 
% Manchester U.K., Springer Editorial, 2nd Ed. Chapert 11.
% 2. Anand Dhariya; 1D Euler Implemetation with shock. Michiga U. 2007;
% http://sitemaker.umich.edu/anand/files/project_2.pdf
%
% This code solves the Sod's shock tube problem 
%
% t=0                                 t=tEnd
% Density                             Density
%   ****************|                 *********\
%                   |                           \
%                   |                            \
%                   |                             ****|
%                   |                                 |
%                   |                                 ****|
%                   ***************                       ***********
%
% coded by Manuel Diaz, 2012.12.25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; %close all; clc;

%% Parameters
CFL   = 0.55;	% CFL number
tEnd  = 0.2;	% Final time
nE    = 100;    % Number of cells/Elements
n     = 5;      % Number of degrees of freedom
alpha = 0.00;	% Parameter for artificial viscosity
IC    = 01;     % 12 IC cases are available

% Ratio of specific heats 
gamma = (n+2)/n;

% Discretize spatial domain
a=0; b=1; dx=(b-a)/nE; nx=nE+1; x=linspace(a,b,nx);

% Set IC
[r0,u0,p0,~,~] = Euler_IC1d(x,IC);
E0 = p0./((gamma-1)*r0)+0.5*u0.^2;  % Total Energy
a0 = sqrt(gamma*p0./r0);            % Speed of sound

% Discretize time domain
dt=CFL*dx/max(abs(u0+a0));  % using the system's largest eigenvalue
t = 0:dt:tEnd; 

% Exact solution
[xe,re,ue,pe,ee,te,Me,se] = ...
   EulerExact(r0(1),u0(1),p0(1),r0(nx),u0(nx),p0(nx),tEnd,n);

%% Solver Loop
% Load initial condition
r=r0; u=u0; p=p0; E=E0; it=0;

visc = zeros(3,nx-1);
for tsteps=t
    % iteration counter
    it=it+1;
    
    % define vectors q & F for every x(i)
    q=[r; r.*u; r.*E];
    F=[r.*u; r.*u.^2+p; u.*(r.*E+p)];
    
    %Calculate q* and flow properties
    q_star(:,1:nx-1) = 0.5*(q(:,1:nx-1)+q(:,2:nx))...
                        - dt/(2*dx)*(F(:,2:nx)-F(:,1:nx-1));
    r(1:nx-1)=q_star(1,1:nx-1);
    u(1:nx-1)=q_star(2,1:nx-1)./r(1:nx-1);
    E(1:nx-1)=q_star(3,1:nx-1)./r(1:nx-1);
    p(1:nx-1)=(gamma-1)*r(1:nx-1).*(E(1:nx-1)-0.5*u(1:nx-1).^2);
    
    %Calculate F*
    F_star(1:3,1:nx-1) = [r(1:nx-1).*u(1:nx-1);...
                        r(1:nx-1).*u(1:nx-1).^2+p(1:nx-1);...
                        u(1:nx-1).*(r(1:nx-1).*E(1:nx-1)+p(1:nx-1))];
    mx = [zeros(1,nx-1);ones(1,nx-1);u(1:nx-1)];
    
    %Calculate arfiticial viscosity
    for j=1:3
        visc(j,1:nx-1) = alpha*dx^2*r(1:nx-1)...
                        .*abs((u(2:nx)-u(1:nx-1))/dx)...
                        .*((u(2:nx)-u(1:nx-1))/dx).*mx(j,1:nx-1);
    end
    F_star=F_star-visc;
    
    %Update F* and q matrix
    q(1:3,2:nx-1)=q(1:3,2:nx-1)-dt/dx*(F_star(1:3,2:nx-1)-F_star(1:3,1:nx-2));
    
    % compute flow properties
    r=q(1,:);  u=q(2,:)./r;  E=q(3,:)./r;  p=(gamma-1)*r.*(E-0.5*u.^2);
end

% Calculation of flow parameters
a = sqrt(gamma*p./r);
M = u./a;
p_ref = 101325;	% Reference air pressure (N/m^2)
rho_ref= 1.225;	% Reference air density (kg/m^3)
s = 1/(gamma-1)*(log(p/p_ref)+gamma*log(rho_ref./r)); 
                % Entropy w.r.t reference condition
ss = log(p./r.^gamma);
                % Dimensionless Entropy
Q = r.*u;     % Mass Flow rate per unit area
e = p./((gamma-1)*r); % internal Energy

%% Final plot
offset=0.05;
s1=subplot(2,3,1); plot(x,r,'.b',xe,re,'k'); xlabel('x(m)'); ylabel('Density (kg/m^3)');
s2=subplot(2,3,2); plot(x,u,'.m',xe,ue,'k'); xlabel('x(m)'); ylabel('Velocity (m/s)');
s3=subplot(2,3,3); plot(x,p,'.k',xe,pe,'k'); xlabel('x(m)'); ylabel('Pressure (Pa)');
s4=subplot(2,3,4); plot(x,ss,'.c',xe,se,'k'); xlabel('x(m)'); ylabel('Entropy/R gas');
s5=subplot(2,3,5); plot(x,M,'.g',xe,Me,'k'); xlabel('x(m)'); ylabel('Mach number');
s6=subplot(2,3,6); plot(x,e,'.r',xe,ee,'k'); xlabel('x(m)'); ylabel('Internal Energy (kg/m^2s)');
title(s1,'Lax-Wendroff Euler Solver');