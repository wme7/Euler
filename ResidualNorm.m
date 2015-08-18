function [norms] = ResidualNorm(res,vol,nN)
%**************************************************************************
% This subroutine computes the residual norms: L1, L2, L_infty
% -------------------------------------------------------------------------
%  Input:  node(:).res = the residuals of every node.
%                  vol = dual volume of every node.
%                   nN = total number of nodes.
% Output:     res_norm = residual norms (L1, L2, Linf).
% -------------------------------------------------------------------------
% NOTE: It is not done here, but I advise you to keep the location of the
%       maximum residual (L_inf).
%**************************************************************************

% Allocate norms arrays
res_norm(:,1) =  0;     % L1   norm
res_norm(:,2) =  0;     % L2   norm
res_norm(:,3) = -1;     % Linf norm

% Compute norms
for i = 1, nnodes
    residual = abs( res./vol ); % Divided residual
    res_norm(:,1) = res_norm(:,1) + residual;       
    res_norm(:,2) = res_norm(:,2) + residual^2;     
    res_norm(:,3) = max(res_norm(:,3), residual);   
end

% Average Residuals Norms
res_norm(:,1) =      res_norm(:,1)/nN;
res_norm(:,2) = sqrt(res_norm(:,2)/nN);

% Norms
norms = res_norm;