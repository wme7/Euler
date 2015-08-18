function gradw = LSQgradients2d(xi,yi,wi,xk,yk,wk,invA)
%**************************************************************************
% Given the primitive variables at k-nodes around node i, this computes the 
% gradients (wx,wy) at node j by the unweighted least-squares method.
%
% -------------------------------------------------------------------------
%  Input:    w  = current primitive variables at node 'inode'
% Output: gradw = [wx,wy] gradients of the primitive variables
% -------------------------------------------------------------------------
%  
%  1. At node i, compute a vector: 
%           b = \sum_k [ (xk-xi)*(wk-wi), (yk-yi)*(wk-wi) ]
%  2. Compute the gradient by multiplying b by the inverse LSQ matrix that
%     has been pre-computed by the subroutine lsq01_matrix() in the main:
%             wx = invA(1,1)*b(1) + invA(1,2)*b(2)
%             wy = invA(2,1)*b(1) + invA(2,2)*b(2)
%
%**************************************************************************

% Compute vector b
b = [ sum((xk-xi).*(wk-wi)), sum((yk-yi).*(wk-wi)) ];
    
% Multiply the inverse LSQ matrix to get the gradients
wx = invA(1,1)*b(1)+invA(1,2)*b(2);
wy = invA(2,1)*b(1)+invA(2,2)*b(2);

% the gradient
gradw = [wx,wy];
    
end