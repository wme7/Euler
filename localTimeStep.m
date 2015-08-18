function dt = localTimeStep(node)
%*******************************************************************************
% This subroutine computes the explicit time-step: the minimum dt over nodes.
%
% ------------------------------------------------------------------------------
%  Input: node(i)%vol = Dual volume
%         node(i)%wsn = Sum of the max wave speed multiplied by the face length
%
% Output:         dt  = global time step
%         node(:)%dt  =  local time step
% ------------------------------------------------------------------------------
%
% NOTE: Local time step is computed and stored at every node, but not used.
%       For steady problems, it can be used to accelerate the convergence.
%
%*******************************************************************************
 subroutine compute_time_step(dt)

 use constants   , only : p2, half, one, two
 use my_main_data, only : nnodes, node

 implicit none

 real(p2), intent(out) ::  dt
%Local variables
 integer  :: i
 real(p2) :: dt_min

 dt_min = 1.0e+05_p2

%--------------------------------------------------------------------------------
  nodes : do i = 1, nnodes
%--------------------------------------------------------------------------------

% Local time step: dt = volume/sum(0.5*max_wave_speed*face_area).

    node(i)%dt = node(i)%vol / node(i)%wsn

% Keep the minimum dt

   if (i==1) dt_min = node(i)%dt
   dt_min = min( dt_min, node(i)%dt )

%--------------------------------------------------------------------------------
  end do nodes
%--------------------------------------------------------------------------------


% Global time-step

   dt = dt_min;
