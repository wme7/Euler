function [x,tRatio] = CouetteFlow(varargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Exact Solution of Navier-Stokes Couette Flow
%                  coded by Manuel Diaz, NTU, 2015.03.21
%
%             (T-T0)/(T1-T0) = (y/H)+ (Pr Ec)/2*(y/H)*(1-y/H)
%
% Refs:
% [1] Xu, Kun. "A gas-kinetic BGK scheme for the Navier-Stokes equations
%     and its connection with artificial dissipation and Godunov method."
%     Journal of Computational Physics 171.1 (2001): 289-335.  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 0
        Ec = 40;
        Pr = 0.5;
        H  = 1.0;
    case 2
        Ec = varargin{1};
        Pr = varargin{2};
        H  = 1.0;
    case 3
        Ec = varargin{1};
        Pr = varargin{2};
        H  = varargin{3};
end
y = 0:0.01:H; x = y/H;
% Compute ration of temperatures: (T-T0)/(T1-T0)
tRatio = x + Pr*Ec/2.*x.*(1-x);