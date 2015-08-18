%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Steady shock structure based on the Navier-Stokes equations.
%               coded by Manuel Diaz, NTU, 14.05.2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Refs:
%   [1] Xu, Kun. "A gas-kinetic BGK scheme for the Navier?Stokes equations
%       and its connection with artificial dissipation and Godunov method."
%       Journal of Computational Physics 171.1 (2001): 289-335. 
%   [2] Masatsuka, K. "I do like CFD." Published by Katate Masatsuka 2009.
%       http://www.cfdbooks.com/ 
% 
% Following implementation in [2] and extended by Manuel Diaz, 14.05.2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Notes: 
%
%    uL(>1) **********     Shock width depends on the Renolds number.
%                     *
%                      *
%                       *********** uR
% 
%  This code solves the 1D Navier-Stokes equations across a steady shock.
%  The steady shock is defined in Section 7.11.8 of "I do like CFD, VOL.1".
% 
%  Input ---------------------------------------------------
%     nnodes = # of nodes on a grid.
%      gamma = Ratio of specific heats (1.4 for air)
%        Pra = Prandtl number
%      M_inf = Upstream Mach number
%     Re_inf = Upstream Reynolds number
%      T_inf = Upstream temperature (K)
%       xmin = Left end of the domain
%       xmax = Right end of the domain
%        eps = Perturbation to the initial velocity.
%  
%  Output --------------------------------------------------
%          x = node coordinate
%        rho = density
%          u = velocity
%          p = pressure
%        tau = viscous stress
%          q = heat flux
% 
%   Note: the solutions are nondimensionalized values (compressible scale),
%           rho=rho/rho_inf, rho=u/a_inf, rho=p/(rho_inf*a_inf^2), 
%             T=T/T_inf    ,  mu=mu/mu_inf.
% 
%   Note: the shock location depends on 'eps'. If eps=0, nothting will
%         happen and the computation of the viscous shock fails. The shock
%         appear sooner for larger value of eps. 
% 
%   Note: The exact solutions are point values. For high-order (>2nd)
%         cell-centered codes (that compute cell-averages, you need to
%         modify this code to output cell-averages.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; %clc; close all;

global gamma Pra C M_inf Re_inf T_inf

%% 0. Input parameters
      nx =  81;     % Number of nodes in the grid
       n =  5;      % Gas total degrees of freedom ( for air)
     Pra =  1;%3/4;	% Prandtl number
    xmin = -1.0;	% Left end of the domain
    xmax =  1.0;	% Right end of the domain
     eps = 1e-12;	% Perturbation to the initial velocity.

%% 1. Generate a grid
dx = (xmax-xmin)/(nx-1); % Domain length divided by ncells = nnodes-1.
 x = linspace(xmin,xmax,nx);

%% 2. Set up the left and right states of the shock (Rankine-Hugoniot).
    gamma = (n+2)/n;	% Ratio of specific heats (1.4 for air)
    M_inf = 1.0;     % Upstream Mach number
    Re_inf= 25.0;	% Upstream Reynolds number
    T_inf = 400.0;	% Upstream temperature (K)
    ML = 1.8;

% Left state
    rL = 1;
    uL = 1;
    %pL = 1/(gamma);
    pL = 1/(gamma*ML^2);
    TL = gamma*pL/rL;
% Right state
    rR = rL*(gamma+1)*ML^2 / ( (gamma-1)*ML^2 + 2 );
    uR = uL/rR;
    %pR = (1+2*gamma/(gamma+1)*(ML^2-1))/(gamma*ML^2);
    pR = (2*gamma/(gamma+1)*ML^2-(gamma-1)/(gamma+1))*pL;
    TR = gamma*pR/rR;
      
% Contansts from the left and right states
% Left state:
   C(1) = rL*uL;                                  	%Mass flux
   C(2) = rL*uL^2 + pL;                             %Momentum flux
   C(3) = uL*( gamma*pL/(gamma-1) + 0.5*rL*uL^2 );	%Energy flux
% Right state:
   C(1) = rR*uR;                                  	%Mass flux
   C(2) = rR*uR^2 + pR;                          	%Momentum flux
   C(3) = uR*( gamma*pR/(gamma-1) + 0.5*rR*uR^2 );	%Energy flux
   % (we can use either state as they are numerically the same!)
% The right state is selected as we integrate from x(N) ---> x(1).

%% 3. Set initial values for the ODE: 
% the right end of the domain (the right end of the grid).
     r(nx) = rR;
     u(nx) = uR;
     p(nx) = pR;
     T(nx) = TR;
   tau(nx) = 0;
     q(nx) = 0;

%% 4. Loop over nodes: 
% compute the exact solution at each node by integrating the ODE for u and
% T from the previous node. 

%Solving the ODE from one node to the next
for i = nx-1:-1:1
    % Previous node location and the solution there.
    x0 = x(i+1);
    u0 = u(i+1);
    T0 = T(i+1);
    
    % The current node location where we compute the solution.
    xi = x(i);

    % Perturbation to the initial velocity. 
    % ( This is required to generate the viscous shock. Without this,
    % nothing happens and you'll never reach the left state ).
    if i==nx-1; u0=u0*(1+eps); end

    % Integrate ODE from x0 to xi.
    [ui,Ti] = ODEsolverRK4(x0,u0,T0,xi);

    % Store the computed solution at the current node (x=xi).
      r(i) = C(1)/ui;
      u(i) = ui;
      p(i) = r(i)*Ti/gamma;
      T(i) = Ti;
    tau(i) = ( C(1)*ui + (C(1)/ui)*(Ti/gamma) - C(2) );                   
      q(i) = ( C(1)*ui^2/2-C(1)/(gamma-1)*(Ti/gamma)+C(3)-C(2)*ui );
end 

%% 5. Save results to a mat file
% Open a Matlab folders to store results
folder = 'CNSshockExact'; t=T;

% save macroscopic variables into their folders
mkdir(folder,'rho');  save([folder,'/rho/rho.mat'],'r');
mkdir(folder,'ux');   save([folder,'/ux/ux.mat'],'u');
mkdir(folder,'p');    save([folder,'/p/p.mat'],'p');
mkdir(folder,'t');    save([folder,'/t/t.mat'],'t');
mkdir(folder,'tau');  save([folder,'/tau/tau.mat'],'tau');
mkdir(folder,'q');    save([folder,'/q/q.mat'],'q');
mkdir(folder,'x');    save([folder,'/x/x.mat'],'x');

%% 6. Visualize solution Profiles
figure(1);
subplot(2,3,1); plot(x,r,'r*'); axis('square'); grid; title('\rho'); xlabel('x'); ylabel('Density \rho'); 
subplot(2,3,2); plot(x,u,'r*'); axis('square'); grid; title('u'); xlabel('x'); ylabel('Velocity u'); 
subplot(2,3,3); plot(x,p,'r*'); axis('square'); grid; title('p'); xlabel('x'); ylabel('Pressure p'); 
subplot(2,3,4); plot(x,t,'r*'); axis('square'); grid; title('T'); xlabel('x'); ylabel('Temperature T'); 
subplot(2,3,5);plot(x,tau,'r*');axis('square'); grid; title('\tau'); xlabel('x'); ylabel('Viscous stress \tau'); 
subplot(2,3,6); plot(x,q,'r*'); axis('square'); grid; title('q'); xlabel('x'); ylabel('Heat flux q'); 