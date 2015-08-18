%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 1-D Compresssible Navier-Stokes Equations
%                      by Manuel Diaz, NTU, 20.03.2015
%
% A finite volume implementaion for the compresssible navier-stokes (CNS)
% equations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% [1] Are Skoinen, Cartesian Grid methods for the compressible
% Navier-Stokes equations. Master Thesis, NTNU Norway, Spring 2012.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; %clc; close all;

%% Numerical Parameters
CFL     = 0.60;     % CFL number
VNN     = 0.40;     % Von Neumann condition
tEnd    = 0.30;     % Final time
nx      = 500;      % Number of cells/Elements in x
n       = 5;        % Number of degrees of freedom
IC      ='Sod';     % 19 IC cases are available
limiter ='MC';      % MM, MC, VA.
fluxMth ='HLLE';    % LF, RUS, ROE, HLLE, HLLC(!).
plot_fig= 0;        % 1:visualize evolution 

%% Physical parameters
% Discretize spatial domain
Lx=1; dx=2*Lx/nx; xc=-Lx+dx/2:dx:Lx;
 gamma = (n+2)/n;   % Ratio of specific heats for ideal di-atomic gas
 M_inf = 1.0;       % Upstream Mach number [-]
Re_inf = 200.0; 	% Upstream Reynolds number [-]
 T_inf = 400.0;     % Upstream temperature [k]
 L_inf = Re_inf;	% Characteristic Length
    Pr = 1;%3/4;	% Prandtl number [-]
    ML = 1.8;

%% Preprocess

% Set IC
[r0,u0,p0] = CNS_IC1d(xc,IC,n,ML);
E0 = p0./((gamma-1)*r0)+0.5*u0.^2;  % Total Energy
c0 = sqrt(gamma*p0./r0);            % Speed of sound
Q0 = [r0; r0.*u0; r0.*E0];
T0 = gamma*p0./r0;

% Compute transport coefs
 C = 110.5; %[K] : for air
mu0= (1+C/T_inf)./(T0+C/T_inf).*T0.^(3/2) * M_inf/Re_inf; 
k0 = gamma*mu0/(Pr*(gamma-1)) * M_inf/Re_inf/gamma;

% Exact solution from invicid theory
% [xe,re,ue,pe,ee,te,Me,se] = ...
%    EulerExact(r0(1),u0(1),p0(1),r0(nx),u0(nx),p0(nx),tEnd,n);

% Set q-array & adjust grid for ghost cells
nx=nx+2; zero=[0;0;0]; q0=[zero,Q0,zero];

% Boundary Conditions in ghost cells
q0(:,1)=q0(:,2); q0(:,nx)=q0(:,nx-1);   % Natural BCs

% Initial time step
lambda0=abs(u0)+c0; a0= max(lambda0(:));
dtInvicid = CFL*dx/a0;
lambda3=4/3*mu0./r0; lambda4=gamma*mu0./(r0*Pr); b0=max([lambda3(:);lambda4(:)]);
dtViscous = VNN/((1/dx^2)*b0);
dt0=min(dtInvicid,dtViscous);

%% Solver Loop

% load initial conditions
q=q0; t=0; it=0; dt=dt0; a=a0; b=b0;

% solve loop
tic
while t < tEnd
    % RK3 step 1
    qs = q - dt*MUSCL_EulerSys(q,a,gamma,dx,nx,limiter,fluxMth) ...
        + dt*fluxCNS1d(q,gamma,dx,nx,Pr,T_inf,Re_inf,M_inf,L_inf);
    
    qs(:,1)=qs(:,2); qs(:,nx)=qs(:,nx-1);   % Natural BCs
    
    % RK3 step 2
    qs2 = (3*q + qs - dt*MUSCL_EulerSys(qs,a,gamma,dx,nx,limiter,fluxMth) ...
        + dt*fluxCNS1d(qs,gamma,dx,nx,Pr,T_inf,Re_inf,M_inf,L_inf))/4;
    
    qs2(:,1)=qs2(:,2); qs2(:,nx)=qs2(:,nx-1);   % Natural BCs
    
    % RK3 step 3 / update solution
    q = (q + 2*(qs2 - dt*MUSCL_EulerSys(qs2,a,gamma,dx,nx,limiter,fluxMth) ...
        + dt*fluxCNS1d(qs2,gamma,dx,nx,Pr,T_inf,Re_inf,M_inf,L_inf)))/3;
    
    q(:,1)=q(:,2); q(:,nx)=q(:,nx-1);   % Natural BCs
    
    % compute flow properties
    r=q(1,:);u=q(2,:)./r;E=q(3,:)./r;p=(gamma-1)*r.*(E-0.5*u.^2);
    c=sqrt(gamma*p./r); T=gamma*p./r;
    
    % compute transport coefs
    mu = transportCoefs(q,gamma,Pr,T_inf,Re_inf,M_inf,L_inf);
    
    % Update dt and time
    lambda0=abs(u)+c; a=max(lambda0(:));
    lambda3=4/3*mu./r; lambda4=gamma*mu./(r*Pr); b=max([lambda3(:);lambda4(:)]);
    dt = min(CFL*dx/a,VNN/((1/dx^2)*b));
    if t+dt>tEnd; dt=tEnd-t; end; t=t+dt; it=it+1;
    
    % Plot figure
    if rem(it,10) == 0
        if plot_fig == 1;
            subplot(2,2,1); plot(xc,r(2:nx-1),'.b');
            subplot(2,2,2); plot(xc,u(2:nx-1),'.m');
            subplot(2,2,3); plot(xc,p(2:nx-1),'.k');
            subplot(2,2,4); plot(xc,T(2:nx-1),'.r');
            drawnow
        end
    end
end
cputime = toc;

% Remove ghost cells
q=q(:,2:nx-1); nx=nx-2; 

%% Compute flow properties
r=q(1,:); u=q(2,:)./r; E=q(3,:)./r; p=(gamma-1)*r.*(E-0.5*u.^2); T=gamma*p./r;

% compute transport coefs
[mu,k] = transportCoefs(q,gamma,Pr,T_inf,Re_inf,M_inf,L_inf);

% heat flux and share stress
tau = 4/3*mu.*gradient(u,dx);
qx = -k.*gradient(T,dx);

%% Save results to a mat file
% Open a Matlab folders to store results
folder = 'CNSshockMUSCL'; t=T; x=xc; q=qx;

% Save macroscopic variables into their folders
mkdir(folder,'rho');  save([folder,'/rho/rho.mat'],'r');
mkdir(folder,'ux');   save([folder,'/ux/ux.mat'],'u');
mkdir(folder,'p');    save([folder,'/p/p.mat'],'p');
mkdir(folder,'t');    save([folder,'/t/t.mat'],'t');
mkdir(folder,'tau');  save([folder,'/tau/tau.mat'],'tau');
mkdir(folder,'q');    save([folder,'/q/q.mat'],'q');
mkdir(folder,'x');    save([folder,'/x/x.mat'],'x');

%% Plots results
figure(2); %title('RK2 TVD-MUSCL Compressible N-S');
subplot(2,3,1); plot(xc,r,'r*'); axis('square'); grid; title('\rho'); xlabel('x'); ylabel('Density \rho'); 
subplot(2,3,2); plot(xc,u,'r*'); axis('square'); grid; title('u'); xlabel('x'); ylabel('Velocity u'); 
subplot(2,3,3); plot(xc,p,'r*'); axis('square'); grid; title('p'); xlabel('x'); ylabel('Pressure p'); 
subplot(2,3,4); plot(xc,T,'r*'); axis('square'); grid; title('T'); xlabel('x'); ylabel('Temperature T'); 
subplot(2,3,5);plot(xc,tau,'r*');axis('square'); grid; title('\tau'); xlabel('x'); ylabel('Viscous stress \tau'); 
subplot(2,3,6); plot(xc,qx,'r*'); axis('square'); grid; title('q'); xlabel('x'); ylabel('Heat flux q'); 