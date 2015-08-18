function [r_0,u_0,v_0,p_0] = CNS_IC2d(x,y,gamma,input)
% Load the IC of a 1D Riemann classical schok tube problem configuration. 
% In the notation we take advantage of the matlab array notation as follows
%
%   1.0 +-----------+-----------+
%       |           |           |       
%       |   reg 2   |   reg 1   |
%       |           |           |
%   0.5 +-----------o-----------+
%       |           |           |
%       |   reg 3   |   reg 4   |
%       |           |           |
%   0.0 +-----------+-----------+
%      0.0         0.5         1.0
%
% prop = [prop_reg1 , prop_reg2 , prop_reg3 , prop_reg4]
%
%   r = rho/density
%   u = velocity in x direction
%   v = velocity in y direction
%   p = Pressure
%
% Manuel Diaz, NTU, 2014.06.27

% Center point (o)
xc = mean(x(:));
yc = mean(y(:));

%% Initial Physical Properties per case:
switch input
    case 'Sod'
        fprintf('Sods Shocktube configuration in 2d \n');
        p = [0.1   1 1 0.1  ];
        r = [0.125 1 1 0.125];
        u = [0     0 0 0    ];
        v = [0     0 0 0    ];
    case 'sShock'
        % Conditions ahead of the shock
        r1 = 1.0;
        u1 = 0.0;
        v1 = 0.0;
        p1 = 0.1;
        % Moving shock relations 
        Ms = 10; Tau = (gamma+1)/(gamma-1); c1 = sqrt(gamma*p1/r1);
        % Conditions behind the shock
        p2 = p1*(2*gamma*Ms^2-(gamma-1))/(gamma+1);
        r2 = r1*(Tau*(p2/p1)+1) / (Tau+(p2/p1));
        u2 = Ms*(1-((gamma-1)*Ms^2+2)/((gamma+1)*Ms^2))*c1;
        v2 = 0.0;
        % Set configuration
        fprintf('Stationary Schock problem configuration in 2d \n');
        p = [p1 p2 p2 p1];
        r = [r1 r2 r2 r1];
        u = [u1 u2 u2 u1];
        v = [v1 v2 v2 v1];
    otherwise
        error('only 18 cases are available');
end
%% Print configuration of selected IC
fprintf('\n');
fprintf('          reg 1 reg 2  reg 3  reg 4\n');
fprintf('density : %2.4f %2.4f %2.4f %2.4f \n',r);
fprintf('  x-vel : %2.4f %2.4f %2.4f %2.4f \n',u);
fprintf('  y-vel : %2.4f %2.4f %2.4f %2.4f \n',v);
fprintf('Presure : %2.4f %2.4f %2.4f %2.4f \n',p);
fprintf('\n');

% Parameters of regions dimensions
reg1 = (x>=xc & y>=yc); % region 1
reg2 = (x <xc & y>=yc); % region 2
reg3 = (x <xc & y <yc); % region 3
reg4 = (x>=xc & y <yc); % region 4

% Initial Condition for our 2D domain
r_0 = r(1)*reg1 + r(2)*reg2 + r(3)*reg3 + r(4)*reg4; % Density, rho
u_0 = u(1)*reg1 + u(2)*reg2 + u(3)*reg3 + u(4)*reg4; % velocity in x
v_0 = v(1)*reg1 + v(2)*reg2 + v(3)*reg3 + v(4)*reg4; % velocity in y
p_0 = p(1)*reg1 + p(2)*reg2 + p(3)*reg3 + p(4)*reg4; % temperature.