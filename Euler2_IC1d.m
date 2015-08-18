function [r0,u0,p0,tEnd,cfl] = Euler2_IC1d(x,input)
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
%
% Refs:
% 1. C.-W. Shu, S. Osher, Efficient implementation of essentially
% non-oscillatory shock-capturing schemes, II, J. Comput. Phys. 83 (1989)
% 32-78. 
% 2. P. Woodward, P.Collela, The numerical simulation of two dimensional
% fluid flow with strong shocks, J. Comput. Phys. 54 (1984) 115-173.
%
%% General Shocktube problems
switch input
    case{1} % Stationary Shock at x=0
        fprintf('Case 1: Stationary Shock, C-W.Shu, S.Osher, JCP 83 (1989) \n');
        rho = [3.857143,0.0,1.0]; % <- modified for sake of this table.
        u   = [-0.920279,-3.549648,-3.549648];
        p   = [10.333,1.0,1.0];
        tEnd = 2.0; cfl = 0.90;
        z1 = (x<=0.0);
        z2 = (x>0.0&x<10);
        z3 = (x>=10);
        % test domain size for this IC
        if (x(1)>-10 || 10>x(end)); error('domain is to small!'); end
        
    case{2} % 
        fprintf('Case 2: Blast Wave, P. Woodward, P.Collela, JCP 53 (1984) \n');
        rho = [1 1 1];
        u   = [0 0 0];
        p   = [1000 0.01 100];
        tEnd = 0.016; cfl = 0.90;
        z1= (x< 0.1);
        z2= (x>=0.1 & x< 0.9);
        z3= (x>=0.9);
        % test domain size for this IC
        if (x(1)>0 || 1>x(end)); error('domain is to small!'); end
        
    otherwise 
        error('Case not available');
        
end
% Print for general Sochtube problems
fprintf('\n');
fprintf('density (z1): %1.3f\n',rho(1));
fprintf('velocity(z1): %1.3f\n',u(1));
fprintf('Presure (z1): %1.3f\n',p(1));
fprintf('\n');
fprintf('density (z2): %1.3f\n',rho(2));
fprintf('velocity(z2): %1.3f\n',u(2));
fprintf('Presure (z2): %1.3f\n',p(2));
fprintf('\n');
fprintf('density (z3): %1.3f\n',rho(3));
fprintf('velocity(z3): %1.3f\n',u(3));
fprintf('Presure (z3): %1.3f\n',p(3));
fprintf('\n');

%% Load Selected case Initial condition:

% Initial Condition for our 2D domain
% Density
if input==1; 
r0 = rho(1).*z1 + (1+0.2*sin(5*x)).*z2 + rho(3).*z3; else
r0 = rho(1).*z1 + rho(2).*z2 + rho(3).*z3; end
% Velovity in x
u0 = u(1).*z1 + u(2).*z2 + u(3).*z3;
% Pressure
p0 = p(1).*z1 + p(2).*z2 + p(3).*z3;

% Plot ICs
%plot(x,r0,'.',x,u0,'.',x,p0,'.'); legend('\rho','u','p',2);