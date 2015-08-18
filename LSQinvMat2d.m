function invAtA = LSQinvMat2d(xi,yi,ni,xk,yk,nk)
%**************************************************************************
% --- Inverse Matrix for 2x2 Least-Squares Gradient Reconstruction ---
%
% Construct a matrix for the linear least-squares(LSQ) gradient reconstruction.
% (unweighted LSQ; more accurate than weighted ones to my knowledge.)
%
% Note: it requires at least 2 non-colinear neighbors.
%
% Example: Consider constructing (ux,uy) at i with the following stencil.
%
%      3 o     o 2
%         \   / 
%          \ /
%         i *-----o 1
%          /|
%         / |
%        /  o 5      *: node in interest (i)
%       o 4          o: neighbors (k = 1,2,3,4,5)
%
%  5 equations:
%    (x1-xi)*ux + (y1-yi)*uy = (u1-ui)
%    (x2-xi)*ux + (y2-yi)*uy = (u2-ui)
%    (x3-xi)*ux + (y3-yi)*uy = (u3-ui)
%    (x4-xi)*ux + (y4-yi)*uy = (u4-ui)
%    (x5-xi)*ux + (y5-yi)*uy = (u5-ui)
%
%  This system is written in the matrix form:
%
%        A*x = b,  x=(ux,uy), A=5x2 matrix, b=5x1 matrix
%
%  The least-squares problem is
%
%      A^T*A*x = A^T*b, (T denotes the transpose: A^T=2x5 matrix)
%  
%  which is
%
%  [sum_k (xk-xi)^2]*ux       + [sum_k (xk-xi)*(yk-yi)]*uy = [sum_k (uk-ui)*(xk-xi)]
%  [sum_k (xk-xi)*(yk-yi)]*ux + [sum_k (yk-yi)]*uy         = [sum_k (uk-ui)*(yk-yi)]
%
% This subroutine computes the inverse of (A^T*A) at every node (which depends
% only on the grid geometry), so that the gradient at a node can be computed
% by a matrix-vector multiplication, i.e., (A^T*A)^{-1}*(A^T*b), 
% (only A^T*b needs to be re-computed).
%
%    if A = | a,b |  then  A^{-1} = 1/(ad-bc) |  d,-b | .
%           | c,d |                           | -c, a |
%
% -------------------------------------------------------------------------
%  Input:    nk  = k neigbour nodes ids
%  Input:    ni  = current node id
%  Input:  xk,yk = neigbour nodes coordinates
%  Input:  xk,yk = current node coordinates
% Output: invAtA = inverse matrix for LSQ reconstruction for current node
% -------------------------------------------------------------------------
%
%**************************************************************************

% Allocate matrix A
a = zeros(2,2);

% Initialize components
a(1,1) = sum((xk - xi).^2);            % dx^2
a(1,2) = sum((xk - xi).*(yk - yi));    % dx*dy
a(2,1) = sum((yk - yi).*(xk - xi));    % dy*dx
a(2,2) = sum((yk - yi).^2);            % dy^2

% Compute determinant
det = a(1,1)*a(2,2) - a(1,2)*a(2,1);
if abs(det)<1E-14; disp([ni,nk]); error('Singular Matrix found!'); end

% Invert and store the inverse matrix:
invAtA = [ a(2,2), -a(2,1); -a(1,2),  a(1,1)]/det;

end
