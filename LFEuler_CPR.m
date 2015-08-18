%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lax-Friedrichs method to solve 1-D Euler equations
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
clear; close all; clc;
global gamma

%% Parameters
CFL     = 0.001;	% CFL number
tEnd    = 0.005;	% Final time
nE      = 50;    % Number of cells/Elements
K       = 1;    % degree of accuaracy
gamma   = 1.4;	% Ratio of specific heats for ideal di-atomic gas
IC      = 8;    % 12 IC cases are available

% Build 1d mesh
xgrid = mesh1d([0 1],nE,'Legendre',K);
dx = xgrid.elementSize; J = xgrid.Jacobian; 
x = xgrid.nodeCoordinates; quad = xgrid.quadratureType;
w = xgrid.weights';	xc = xgrid.elementCenter;

% compute gR'(xi) & gL'(xi)
RR = CorrectionPolynomial('RadauRight',K+1); % g: one-order higher
dg.RR = RR.eval_dP(xgrid.solutionPoints); dg.RL = -flipud(dg.RR);

% Build Lagrange k-Polynomials
l = LagrangePolynomial(xgrid.solutionPoints);
L.lcoef = double(subs(l.lagrangePolynomial,-1));
L.rcoef = double(subs(l.lagrangePolynomial,1));
L.dcoef = double(subs(l.dlagrangePolynomial,xgrid.solutionPoints));

% Set IC
[rho0,u0,p0,~,~] = Euler_IC1d(x,IC);
E0 = p0./((gamma-1)*rho0)+0.5*u0.^2;  % Total Energy
a0 = sqrt(gamma*p0./rho0);            % Speed of sound

% Discretize time domain
dt=CFL*dx/max(max(abs(u0+a0)));  % using the system's largest eigenvalue
t = 0:dt:tEnd; 

% Exact solution
[xe,rhoe,ue,pe,Me,se,ee] = ...
   EulerExact(rho0(1),u0(1),p0(1),rho0(end),u0(end),p0(end),tEnd);

%% Solver Loop

% Low storage Runge-Kutta coefficients
rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0];

% Load initial condition
rho=rho0; u=u0; p=p0; E=E0; it=0;

% Using a 4rd Order 5-stage SSPRK time integration
res_q = zeros(3*(K+1),nE); % Runge-Kutta residual storage

for tsteps=t
    % iteration counter
    it=it+1;
    
    % define vector q for every x in [a,b]
    q=[rho; rho.*u; rho.*E];
    
    % RK's stages
    for RKs = 1:5
        t_local = t + rk4c(RKs)*dt;
        dF = residualCPR(q,K,L,dg);
        res_q = rk4a(RKs)*res_q + dt*dF/J;
        q = q - rk4b(RKs)*res_q;
    end

    % compute flow properties
    rho=q(1:K+1,1:nE);
    u=q(K+2:2*(K+1),1:nE)./rho;
    E=q(2*K+3:3*(K+1),1:nE)./rho;
    p=(gamma-1)*rho.*(E-0.5*u.^2);
%end

%% Calculation of flow parameters
a = sqrt(gamma*p./rho);
M = u./a;
p_ref = 101325;	% Reference air pressure (N/m^2)
rho_ref= 1.225;	% Reference air density (kg/m^3)
s = 1/(gamma-1)*(log(p/p_ref)+gamma*log(rho_ref./rho)); 
                % Entropy w.r.t reference condition
ss = log(p./rho.^gamma);
                % Dimensionless Entropy
Q = rho.*u;     % Mass Flow rate per unit area
e = p./((gamma-1)*rho); % internal Energy

%% Final plot
offset=0.05;
s1=subplot(2,3,1); plot(x,rho,'-+',xe,rhoe,'k'); xlabel('x(m)'); ylabel('Density (kg/m^3)');
s2=subplot(2,3,2); plot(x,u,'-+',xe,ue,'k'); xlabel('x(m)'); ylabel('Velocity (m/s)');
s3=subplot(2,3,3); plot(x,p,'-+',xe,pe,'k'); xlabel('x(m)'); ylabel('Pressure (Pa)');
s4=subplot(2,3,4); plot(x,ss,'-+',xe,se,'k'); xlabel('x(m)'); ylabel('Entropy/R gas');
s5=subplot(2,3,5); plot(x,M,'-+',xe,Me,'k'); xlabel('x(m)'); ylabel('Mach number');
s6=subplot(2,3,6); plot(x,e,'-+',xe,ee,'k'); xlabel('x(m)'); ylabel('Internal Energy (kg/m^2s)');
title(s1,'Lax-Friedrichs Euler Solver');
    
drawnow
end