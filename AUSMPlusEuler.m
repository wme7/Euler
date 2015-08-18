%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUSM+ method to solve 1-D Euler equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Following the ideas of:
% 1. Meng-Sing Liou and Christopher J. Steffen Jr. A New Flux Splitting
%     Scheme. JCP, vol. 10, 23-39 (1991). 
% 2. Meng-Sing Liou. A Sequel to AUSM: AUSM+. JCP, vol. 129, pp 364-382
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
clear all; %close all; clc;
global gamma

%% Parameters
CFL     = 0.5;	% CFL number
tEnd    = 0.14;	% Final time
nE      = 200;  % Number of cells/Elements
gamma   = 1.4;	% Ratio of specific heats for ideal di-atomic gas
IC      = 1;    % 12 IC cases are available
plot_fig= 1;

% Discretize spatial domain
a=0; b=1; x=linspace(a,b,nE); dx=(b-a)/nE;

% Set IC
[rho0,u0,p0,~,~] = Euler_IC1d(x,IC);
E0 = p0./((gamma-1)*rho0)+0.5*u0.^2;  % Specific Total Energy
a0 = sqrt(gamma*p0./rho0);            % Speed of sound

% Discretize time domain
dt=CFL*dx/max(abs(u0+a0));  % using the system's largest eigenvalue
t = 0:dt:tEnd; 

% Exact solution
% [xe,rhoe,ue,pe,Me,se,ee] = ...
%    EulerExact(rho0(1),u0(1),p0(1),rho0(nE),u0(nE),p0(nE),tEnd);

%% Solver Loop
% Load initial condition
rho = rho0; u = u0; p = p0; E=E0; it=0; q_next = zeros(3,nE);
beta = 1/8;  alpha = 3/16;

% M = -1.5:0.1:1.5; 
Mp = @(b,m) 0.*(m<=-1) + ( 1/4*((m+1).^2) +b*(m.^2-1).^2).*(m>-1 & m<1) + m.*(m>=1);
Mm = @(b,m) m.*(m<=-1) + (-1/4*((m-1).^2) -b*(m.^2-1).^2).*(m>-1 & m<1) + 0.*(m>=1);
	
Pp = @(a,m) 0*(m<=-1) + (1/4*(m+1).^2.*(2-m) +a*m.*(m.^2-1).^2).*(m>-1 & m<1) + 1*(m>=1);
Pm = @(a,m) 1*(m<=-1) + (1/4*(m-1).^2.*(2+m) -a*m.*(m.^2-1).^2).*(m>-1 & m<1) + 0*(m>=1);

% For testing: visualize correct behavior of splitting functions
% plot(M,Mp(beta,M),M,Mm(beta,M));
% plot(M,Pm(beta,M),M,Pp(beta,M));

for tsteps=t
    % iteration counter
    it=it+1;
    
    % Define: Left and right masks
    L = 1:nE-1;     R = 2:nE;
    
    % define vectors q & F for every x(i)
    q=[rho; rho.*u; rho.*E];
    F=[rho.*u; rho.*u.^2+p; u.*(rho.*E+p)]; 
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
      
	% Compute M_{j+1/2}^{+} and M_{j+1/2}^{-}
    Mjhalf = Mp(beta,ML) + Mm(beta,MR);
    Mpjhalf = 0.5*(Mjhalf + abs(Mjhalf));
    Mmjhalf = 0.5*(Mjhalf - abs(Mjhalf));
    
    % Compute p_{j+1/2}
    Pjhalf = Pp(alpha,ML).*pL + Pm(alpha,ML).*pR;
    
    % Compute a_{j+1/2}
    ajhalf = 0.5*(aL + aR);
    
    %Positive Part of Flux evaluated in the left cell.
     Fp(1,:) = ajhalf.*(Mpjhalf.*rhoL     + Mmjhalf.*rhoR    );
     Fp(2,:) = ajhalf.*(Mpjhalf.*rhoL.*uL + Mmjhalf.*rhoR.*uR) + Pjhalf;
     Fp(3,:) = ajhalf.*(Mpjhalf.*rhoL.*HL + Mmjhalf.*rhoR.*HR);
   
    % Compute next time step
    for i = 2:nE-1
        q_next(:,i) = q(:,i) - dt/dx*(Fp(:,i) - Fp(:,i-1));
    end
    
    % BCs
    q_next(:,1)  = q(:,2); % Neumann condition to the left
    q_next(:,nE) = q(:,nE-1); % Neumann condition to the right

    % compute flow properties
    rho=q_next(1,:);
    u=q_next(2,:)./q_next(1,:);
    E=q_next(3,:)./q_next(1,:);
    p=(gamma-1)*rho.*(E-0.5*u.^2);
    a = sqrt(gamma*p./rho);
    M = u./a;
    e = E-1/2*u.^2;
    
	% Plot figure
    if plot_fig == 1;
        s1=subplot(2,3,1); plot(x,rho,'.b'); xlabel('x(m)'); ylabel('Density, $\rho(x,t)$','Interpreter','latex')
        s2=subplot(2,3,2); plot(x,u,'.m');  xlabel('x(m)'); ylabel('Velocity, $u(x,t)$','Interpreter','latex')
        s3=subplot(2,3,3); plot(x,p,'.k');  xlabel('x(m)'); ylabel('Pressure, $p(x,t)$','Interpreter','latex')
        s4=subplot(2,3,4); plot(x,M,'.k');  xlabel('x(m)'); ylabel('Mach, $M(x,t)$','Interpreter','latex')
        s5=subplot(2,3,5); plot(x,e,'.b');  xlabel('x(m)'); ylabel('specific Internal energy, $E(x,t)$','Interpreter','latex')
        s6=subplot(2,3,6); plot(x,E,'.r');  xlabel('x(m)'); ylabel('Total specific energy, $E(x,t)$','Interpreter','latex')
        title(s1,'AUSM Euler Solver');
    end
        
    %if rem(it,10) == 0
        drawnow
    %end
end