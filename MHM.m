%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lax-Friedrichs method to solve 1-D Euler equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Following the ideas of:
% 1. E.F. Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics 
%    Manchester U.K., Springer Editorial, 2nd Ed. Chapert 11.
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
clear; close all; clc;

%% Parameters
CFL	= 0.3;	% CFL number
tEnd= 0.1;	% Final time
nE	= 200;  % Number of cells/Elements
n	= 5;	% Number of degrees of freedom
IC	= 08;	% 12 IC cases are available
draw= 1;

% Ratio of specific heats for ideal di-atomic gas
gamma=(n+2)/n;

% Discretize spatial domain
a=0; b=1; x=linspace(a,b,nE); dx=(b-a)/nE;

% Set IC
[rho0,u0,p0,~,~] = Euler_IC1d(x,IC);
E0 = p0./((gamma-1)*rho0)+0.5*u0.^2; % Total Energy
a0 = sqrt(gamma*p0./rho0);           % Speed of sound
q0 = [rho0; rho0.*u0; rho0.*E0];     % Initial condition

% Discretize time domain
dt0=CFL*dx/max(abs(u0+a0));  % using the system's largest eigenvalue

% Exact solution
[xe,rhoe,ue,pe,ee,te,Me,se] = ...
   EulerExact(rho0(1),u0(1),p0(1),rho0(nE),u0(nE),p0(nE),tEnd,n);

%% Solver Loop

% Allocate flux array
HLLE=zeros(3,nE);

% Load initial condition
it=0; t=0; dt=dt0; q=q0;

while t<tEnd
    % iteration counter
    it=it+1; t=t+dt;
        
    % update q matrix and flow parameters
    for i=1:nE-1; HLLE(:,i)=fluxHLLE1d(q(:,i),q(:,i+1),gamma); end
    q(:,2:nE-1)=q(:,2:nE-1) - dt/dx*(HLLE(:,2:nE-1)-HLLE(:,1:nE-2));    
    
	% compute flow properties
    rho=q(1,:);
    u=q(2,:)./rho;
    E=q(3,:)./rho;
    p=(gamma-1)*rho.*(E-0.5*u.^2);
    
    % Plot figure
    if rem(it,10) == 0
        if draw == 1;
            subplot(2,2,1); plot(x,rho,'.b');
            subplot(2,2,2); plot(x,u,'.m'); 
            subplot(2,2,3); plot(x,p,'.k'); 
            subplot(2,2,4); plot(x,E,'.r');
        end
	drawnow
    end
end

% Calculation of flow parameters
a = sqrt(gamma*p./rho); M = u./a;
p_ref = 101325;         % Reference air pressure (N/m^2)
rho_ref= 1.225;         % Reference air density (kg/m^3)
s = 1/(gamma-1)*(log(p/p_ref)+gamma*log(rho_ref./rho)); 
                        % Entropy w.r.t reference condition
ss= log(p./rho.^gamma); % Dimensionless Entropy
Q = rho.*u;             % Mass Flow rate per unit area
e = p./((gamma-1)*rho); % internal Energy

%% Final plot
offset=0.05;
s1=subplot(2,3,1); plot(x,rho,'.b',xe,rhoe,'k'); xlabel('x(m)'); ylabel('Density (kg/m^3)');
s2=subplot(2,3,2); plot(x,u,'.m',xe,ue,'k'); xlabel('x(m)'); ylabel('Velocity (m/s)');
s3=subplot(2,3,3); plot(x,p,'.k',xe,pe,'k'); xlabel('x(m)'); ylabel('Pressure (Pa)');
s4=subplot(2,3,4); plot(x,ss,'.c',xe,se,'k'); xlabel('x(m)'); ylabel('Entropy/R gas');
s5=subplot(2,3,5); plot(x,M,'.g',xe,Me,'k'); xlabel('x(m)'); ylabel('Mach number');
s6=subplot(2,3,6); plot(x,e,'.r',xe,ee,'k'); xlabel('x(m)'); ylabel('Internal Energy (kg/m^2s)');
title(s1,'MHM Euler Solver');