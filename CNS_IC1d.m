function [r0,u0,p0,tEnd,cfl] = CNS_IC1d(x,input,n,ML)
% Load the IC of a classical 1D Riemann Problems.
%
% By Manuel Diaz 2012.10.24.
% In the notation we take advantage of the matlab array notation as follows
%
% prop = [prop_left , prop_right]
% 
% Notation:
% u   = Velocity in x direction
% p   = Pressure
% rho = Density
% r   = Fugacity
% E   = Enerty
% t   = temperature
%
%% Riemann Problems
switch input
    case 'Sod' % Configuration 1, Sod's Problem
        fprintf('Case 1: Sods problem \n');
        p   = [1    0.1  ];
        u   = [0    0    ];
        rho = [1    0.125];
        tEnd = 0.1; cfl = 0.90;
        
    case 'Lax' % Lax test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997
        fprintf('Lax test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997');
        p   = [3.528 0.571];
        u   = [0.698 0    ];
        rho = [0.445 0.5  ];
        tEnd = 0.15; cfl = 0.90; 
        
    case 'sShock' % Stationary Shock
        gamma = (n+2)/n;
        % Left state
        rL = 1;
        uL = 1;
        %pL = 1/(gamma);
        pL = 1/(gamma*ML^2);
        % Right state
        rR = rL*(gamma+1)*ML^2 / ( (gamma-1)*ML^2 + 2 );
        uR = uL/rR;
        %pR = (1+2*gamma/(gamma+1)*(ML^2-1))/(gamma*ML^2);
        pR = (2*gamma/(gamma+1)*ML^2-(gamma-1)/(gamma+1))*pL;

        % set IC
        fprintf('Stationary Shock');
        p   = [ pL  pR ];
        u   = [ uL  uR ];
        rho = [ rL  rR ];
        tEnd = 1.0; cfl = 0.4;
           
    case 'sShock2' % Stationary Shock
        fprintf('Stationary Shock');
        p   = [ 0.1  0.676 ];
        u   = [ 1.2  0.723966942148760 ];
        rho = [ 1.0  1.657534246575342 ];
        tEnd = 0.1; cfl = 0.25;
        
    otherwise 
        error('Case not available');
        
end
% Print for Riemann Problems
fprintf('\n');
fprintf('density (L): %1.3f\n',rho(1));
fprintf('velocity(L): %1.3f\n',u(1));
fprintf('Presure (L): %1.3f\n',p(1));
fprintf('\n');
fprintf('density (R): %1.3f\n',rho(2));
fprintf('velocity(R): %1.3f\n',u(2));
fprintf('Presure (R): %1.3f\n',p(2));
fprintf('\n');

%% Load Selected case Initial condition:
% Pre-Allocate variables
r0 = zeros(size(x)); 
u0 = zeros(size(x)); 
p0 = zeros(size(x));

% Parameters of regions dimensions
x_middle = (x(end)+x(1))/2;
L = find(x<=x_middle);
R = find(x>x_middle);

% Initial Condition for our 2D domain
% Fugacity
r0(L) = rho(1); % region 1
r0(R) = rho(2); % region 2
% Velovity in x
u0(L) = u(1); % region 1
u0(R) = u(2); % region 2
% temperature
p0(L) = p(1); % region 1
p0(R) = p(2); % region 2