%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Rankine-Hugoniot Conditions for Moving and Stationary Shocks
%                      by Manuel Diaz, NTU, 29.04.2015
%
% Also referred to as Rankine-Hugoniot jump conditions or Rankine-Hugoniot
% relations, describe the relationship between the states on both sides of
% a shock wave in a one-dimensional flow in fluids or a one-dimensional
% deformation in solids.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Refs:
%   [1] Rankine-Hugoniot conditions, available online at:
%       http://en.wikipedia.org/wiki/Rankine-Hugoniot_conditions
%   [2] Athena2D in Fortran 95, available online at:
%       http://www.astro.virginia.edu/VITA/ATHENA/dmr.html
%   [3] J.-Y. Yang, J.-C. Huang, L. Tsuei. Numerical solutions of the
%       nonlinear model Boltzmann Equations, Proc. R. Soc. Lond. A, 448,
%       55-80, (1995). 
%   [4] M. Bennounne, M. Lemou, L. Mieussens, Uniformly stable numerical
%       schemes for the Botlzmann equations preserving the compressible
%       Navier-Stokes asymptotics. J. Compp. Phys. 227, 3781-3803, (2008).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% EXAMPLE 1: Double Mach Reflection Initial Setup conditions
gamma = 1.4;

% Conditions ahead of the shock
r1 = 1.4;
%u1 = 0.0; % <--- notice that it is not required in the computation! :O
p1 = 1.0;

%% Moving shock relations 
Ms = 10; Tau = (gamma+1)/(gamma-1); c1 = sqrt(gamma*p1/r1);

% Conditions behind the shock
p2 = p1*(2*gamma*Ms^2-(gamma-1))/(gamma+1);
r2 = r1*(Tau*(p2/p1)+1) / (Tau+(p2/p1));
u2 = Ms*(1-((gamma-1)*Ms^2+2)/((gamma+1)*Ms^2))*c1;

% the components of the reflected shock
u2x = u2*cosd(30);
u2y =-u2*sind(30);

% displays solution
disp([r2,u2x,u2y,p2]);

% Exact solution of double mach reflection problem, see ref. [2],
r2e = 8.0; u2xe = 8.25*cosd(30); v2ye = -8.25*sind(30); p2e = 116.5; 
disp([r2e,u2xe,v2ye,p2e]);

% are they the same?
disp(abs([r2,u2x,u2x,p2]-[r2e,u2xe,v2ye,p2e])<1E6);

%% EXAMPLE 1: Stationary Shock for the Micro-Macro Decomposition
n=1; gamma = (n+2)/n; % gamma = 3, assumed in Mieussens paper!

% Conditions behind the shock
r1 = 1.0;
u1 = 1.2; % <--- for the stationary shock now is it required! :O
t1 = 0.1;
p1 = r1*t1;

%% Moving shock relations 
%Ms = 2.1748; Tau = (gamma+1)/(gamma-1); c1 = sqrt(gamma*t1);
Ms = 2.2; Tau = (gamma+1)/(gamma-1); c1 = sqrt(gamma*t1);

% Conditions ahead of the stationary shock
p2 = p1*(2*gamma*Ms^2-(gamma-1))/(gamma+1);
r2 = r1*(Tau*(p2/p1)+1) / (Tau+(p2/p1));
% u2 = Ms*(1-((gamma-1)*Ms^2+2)/((gamma+1)*Ms^2))*c1;
u2 = u1/(r2/r1); % from eq.(21) in ref.[1], r2/r1 = u1/u2;

% displays solution
disp([r2,u2,p2/r2]); disp(p2);

% Exact solution of double mach reflection problem
r2e = 1.65; u2e = 0.72; t2e = 0.4; 

% Exact solutions presented in ref. [4],
disp([r2e,u2e,t2e]); 
