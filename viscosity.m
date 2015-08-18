function mu = viscosity(varargin)
% Sutherland's law

% Parameter to approximate air for a wide range of temperatures, 
% between 170 [k] and 1900 [K]. 

% 1. Verify Temperature > 0
T = varargin{1};
if T < 0; error('Negative temperature entered to viscosity()!'); end

% 2. Compute viscosity
switch nargin
    case 1 % normalization with respect to free stream
        % Parameters
        C = 110.5; %[K]
        MUo = 1.716E-5; %[Pa.s]
        To = 273.1;     %[K] : for air
        % Sutherland's law (nondimensional form)
        mu = MUo*(To+C)/(T+C)*(T/To)^(3/2);
    case 2  % normalization with respect to reference value
        % Parameters
        C = 110.5;      %[K] : for air
        T_inf = varargin{2};
        % T_inf = 293.15 [K] : Air at sea-level condition (e.g.) 
        % Sutherland's law (nondimensional form)    
        mu = (1+C/T_inf)/(T+C/T_inf)*T^(3/2);
end