%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Generate 1d Exact solutions for Euler system of equations given an IC.
%
%              coded by Manuel Diaz, NTU, 2014.02.08
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
global gamma; gamma = 5/3;

IC    = 1;                      % IC:1~12. See Euler_IC1d.m
tstep = [0.02 0.04 0.06 0.08 0.1];  % End time - Parameter part of ICs
set   = 'auto' ;               % set number of x points of the solution

% Set Domain,
x = linspace(0,1,100);

% IC for physical space (Macroscopic variables)
[rho0,u0,p0,~,~] = Euler_IC1d(x,IC); % 

% Set name for the selected initial condition
IDname = ['ExactSolution_IC',num2str(IC),'.plt'];

% Open a Files to store the Results
file = fopen(IDname,'w');
    % 'file' gets the handel for the file "case.plt".
    % 'w' specifies that it will be written.
    % similarly 'r' is for reading and 'a' for appending.
fprintf(file, 'TITLE = "%s"\n',IDname);
fprintf(file, 'VARIABLES = "x" "density" "velocity" "temperature" "pressure" "energy" "entropy"\n');

% for every tstep in the list,
for i = 1:length(tstep)
    
    % compute exact solution
    switch set
        case 'auto' % use default size(x)=500
            [xe,re,ue,pe,Me,Se,ee] = ...
                EulerExact(rho0(1),u0(1),p0(1),rho0(end),u0(end),p0(end),tstep(i));
        case 'manual' % set size(x)
            [xe,re,ue,pe,Me,Se,ee] = ...
                EulerExact_modif(rho0(1),u0(1),p0(1),rho0(end),u0(end),p0(end),tstep(i),x);
    end
    te=2/3*ee; nxe=length(xe);
    
    % Write Results
    fprintf(file, 'ZONE T = "time %0.4f"\n', tstep(i));
    fprintf(file, 'I = %d, J = 1, K = 1, F = POINT\n\n',nxe);
    for j = 1:nxe % to the length ot xe
        fprintf(file, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n', ...
            xe(j),re(j),ue(j),te(j),pe(j),ee(j),Se(j));
    end
end

%% Close file with Results
fprintf('Simulation has been completed succesfully!\n')
fclose(file); fprintf('All Results have been saved!\n')