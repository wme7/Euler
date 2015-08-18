function [w_0,w_inf,M_inf] = ShockDiffractionIC()
%  Shock Diffraction Problem:
%
%                             Wall
%                     --------------------
% Post-shock (inflow) |                  |
% (rho,u,v,p)_inf     |->Shock (M_shock) |            o: Corner node
%    M_inf            |                  |
%              .......o  Pre-shock       |Outflow
%                Wall |  (rho0,u0,v0,p0) |
%                     |                  |
%                     |                  |
%                     --------------------
%                           Outflow
%
gamma = 1.4; % assumed for sea-level air.

% Downstream (ahead of shock) conditions
r0 = 1;
u0 = 0;
v0 = 0;
p0 = 1/gamma;

% Shock Speed
M_shock = 5.09;
u_shock = M_shock * sqrt(gamma*p0/r0);

% Upstream (behind the shock) conditions
r_inf = r0 * (gamma+1)*M_shock^2/((gamma-1)*M_shock^2+2);
p_inf = p0 * ( 2*gamma*M_shock^2 - (gamma-1) )/(gamma+1);
u_inf = (1 - r0/r_inf)*u_shock;
M_inf = u_inf / sqrt(gamma*p_inf/r_inf);
v_inf = 0;

% setting the initial condition in every element of the domain will be made
% outside this function.
w_0 = [r0;u0;v0;p0];  w_inf = [r_inf;u_inf;v_inf;p_inf];