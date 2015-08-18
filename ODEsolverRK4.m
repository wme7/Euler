function [ui,Ti] = ODEsolverRK4(x0,u0,T0,xi)
% This subroutine integrates the ODE from x=x0 to x=xi and find the
% velocity and temperature at x=xi: ui and Ti.
%
% Input ---------------------------------------------------
%   x0 =    initial position
%   u0 =    velocity at x=x0
%   T0 = temperature at x=x0
%
% Output --------------------------------------------------
%   xi =      final position
%   ui =    velocity at x=xi
%   Ti = temperature at x=xi
%
% Note: the ODE is integrated from x=x0 to x=xi in 1000 steps by the
%       classical 4th-order Runge-Kutta method.
%
% Increment for ODE integration: 1000 steps between x0 and xi.
dz = (x0-xi)/1000;

% Integrate ODE, dV/dx=RHS, from x0 to xi by the classical RK4.
%      xi <--- x0
%  -----o-------o-------->x

% We integrate actually dV/dz=-RHS where z=-x.
%      z0 ---> zi
%  -----o-------o-------->z

% 1. Initial values and final location, zi, in the reversed coordinate.
z0  = -x0;
V(1)=  u0;
V(2)=  T0;
zi  = -xi;

% 2. Stepping from -x0 to -xi. Solve the ODE with x=-x.
z = z0;

% do z0 to zi:
finish = false;
while ~finish
    % To finish up at z=zi precisely
    if (z + dz > zi); dz = zi - z; finish = true; end
    z = z + dz;
	K1 = V + 0.5*dz*( rhs( V) );
	K2 = V + 0.5*dz*( rhs(K1) );
	K3 = V +     dz*( rhs(K2) );
	V = (K1 + 2*K2 + K3 - V)/3 + dz*( rhs(K3) )/6;
end

% Solution at z=zi, i.e., x=xi.
ui =  V(1);
Ti =  V(2);

end
     
function RHS = rhs(V)

global gamma Pra C M_inf Re_inf T_inf

u=V(1); T=V(2); mu=viscosity(T,T_inf);

% Coefs * Nondimensionalization factor (compressible scale) 
Cu=  4/3*mu*M_inf/Re_inf;
CT= -gamma*mu/(Pra*(gamma-1))*M_inf/Re_inf/gamma;

RHS(1)=-( C(1)*u + C(1)*T/gamma/u - C(2) )/Cu;
RHS(2)=-( C(1)*u^2/2-C(1)/(gamma-1)*T/gamma+C(3)-C(2)*u )/CT;
  
end