%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUSM method to solve 1-D Euler equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Following the ideas of:
% 1. Meng-Sing Liou and Christopher J. Steffen Jr. A New Flux Splitting
%     Scheme. JCP, vol. 10, 23-39 (1991). 
% 2. Meng-Sing Liou. A Sequel to AUSM: AUSM. JCP, vol. 129, pp 364-382
%     (1996). 
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
% coded by Manuel Diaz, 2014.05.14.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; %close all; clc;

%% Parameters
CFL     = 0.5;	% CFL number
tEnd    = 0.2;	% Final time
nE      = 200;  % Number of cells/Elements
n       = 5;	% ideal gas total degrees of freedom
IC      = 1;    % 12 IC cases are available
plot_fig= 1;

% Ratio of specific heats for ideal di-atomic gas
gamma = (n+2)/n;

% Discretize spatial domain
a=0; b=1; x=linspace(a,b,nE); dx=(b-a)/nE;

% Set IC
[r0,u0,p0] = Euler_IC1d(x,IC);
E0 = p0./((gamma-1)*r0)+0.5*u0.^2;  % Specific Total Energy
a0 = sqrt(gamma*p0./r0);            % Speed of sound

% Discretize time domain
dt=CFL*dx/max(abs(u0+a0));  % using the system's largest eigenvalue
t = 0:dt:tEnd; 

% Exact solution
[xe,re,ue,pe,Me,se,ee] = ...
   EulerExact(r0(1),u0(1),p0(1),r0(nE),u0(nE),p0(nE),tEnd,n);

%% Solver Loop
% Load initial condition
r=r0; u=u0; p=p0; E=E0; it=0; q_next=zeros(3,nE);

for tsteps=t
    % iteration counter
    it=it+1;
    
    % Define: Left and right masks
    L = 1:nE-1;     R = 2:nE;
    
    % define vectors q & F for every x(i)
    q=[r; r.*u; r.*E];
    F=[r.*u; r.*u.^2+p; u.*(r.*E+p)]; 
    % this F is not used here, just wrote it for the sake of completeness.
    
    % Left and right primitive variable arrays
    qL = q(:,L);    qR = q(:,R);
    
    % Left state
    rhoL = qL(1,:);
      uL = qL(2,:)./qL(1,:);
      pL = (gamma-1)*(qL(3,:)-1/2*rhoL.*uL.*uL);
      aL = sqrt(gamma*pL./rhoL);
      ML = uL./aL;
      HL = (qL(3,:)+pL)./rhoL;
    % Right state
    rhoR = qR(1,:);
      uR = qR(2,:)./qR(1,:);
      pR = (gamma-1)*(qR(3,:)-1/2*rhoR.*uR.*uR);
      aR = sqrt(gamma*pR./rhoR);
      MR = uR./aR;
      HR = (qR(3,:)+pR)./rhoR;
      
	% Compute Mp and Pp @ Left cell
    %ML = -1.5:0.1:1.5; pL = ones(size(ML)); % for testing
    Mp = 0*(ML<=-1) +...
        (1/4*(ML+1).^2).*(ML>-1 & ML<1) +...
        ML.*(ML>=1);
    Pp = 0*(ML<=-1) +...
        (1/4*pL.*(1+ML).^2.*(2-ML)).*(ML>-1 & ML<1) +...
        pL.*(ML>=1);
    
    % Compute Mm and Pm @ Right cell
    % MR = -1.5:0.1:1.5; pR = ones(size(ML)); % for testing
    Mm = MR.*(MR<=-1) +...
        (-1/4*(MR-1).^2).*(MR>-1 & MR<1) +...
        0*(MR>=1);
    Pm = pR.*(MR<=-1) +...
        (1/4*pR.*(1-MR).^2.*(2+MR)).*(MR>-1 & MR<1) +...
        0*(MR>=1);

    % For testing: visualize correct behavior of splitting functions
    %plot(ML,Mp,ML,Pp,MR,Mm,MR,Pm); 
    
    %Positive Part of Flux evaluated in the left cell.
     Fp(1,:) = max(0,Mp+Mm).*aL.*rhoL;
     Fp(2,:) = max(0,Mp+Mm).*aL.*rhoL.*uL + Pp;
     Fp(3,:) = max(0,Mp+Mm).*aL.*rhoL.*HL;

    %Negative Part of Flux evaluated in the right cell.
     Fm(1,:) = min(0,Mp+Mm).*aR.*rhoR;
     Fm(2,:) = min(0,Mp+Mm).*aR.*rhoR.*uR + Pm;
     Fm(3,:) = min(0,Mp+Mm).*aR.*rhoR.*HR;

    %Compute the flux: Fp(uL)+Fm(uR).
        AUSM = Fp + Fm;
    
    % Compute next time step
    for i = 2:nE-1
        q_next(:,i) = q(:,i) - dt/dx*(AUSM(:,i) - AUSM(:,i-1));
    end
    
    % BCs
    q_next(:,1)  = q(:,2);      % Neumann condition to the left
    q_next(:,nE) = q(:,nE-1);   % Neumann condition to the right

    % compute flow properties
    r=q_next(1,:);
    u=q_next(2,:)./q_next(1,:);
    E=q_next(3,:)./q_next(1,:);
    e = E-0.5*u.^2;
    p=(gamma-1)*r.*e;
    a = sqrt(gamma*p./r);
    M = u./a;
    
	% Plot figure
    if plot_fig == 1;
        s1=subplot(2,3,1); plot(xe,re,x,r,'.b'); xlabel('x(m)'); ylabel('Density, $\rho(x,t)$','Interpreter','latex')
        s2=subplot(2,3,2); plot(xe,ue,x,u,'.m');  xlabel('x(m)'); ylabel('Velocity, $u(x,t)$','Interpreter','latex')
        s3=subplot(2,3,3); plot(xe,pe,x,p,'.k');  xlabel('x(m)'); ylabel('Pressure, $p(x,t)$','Interpreter','latex')
        s4=subplot(2,3,4); plot(x,M,'.k');  xlabel('x(m)'); ylabel('Mach, $M(x,t)$','Interpreter','latex')
        s5=subplot(2,3,5); plot(x,e,'.b');  xlabel('x(m)'); ylabel('specific Internal energy, $E(x,t)$','Interpreter','latex')
        s6=subplot(2,3,6); plot(x,E,'.r');  xlabel('x(m)'); ylabel('Total specific energy, $E(x,t)$','Interpreter','latex')
        title(s1,'AUSM Euler Solver');
    end
        
    if rem(it,10) == 0
        drawnow
    end
end