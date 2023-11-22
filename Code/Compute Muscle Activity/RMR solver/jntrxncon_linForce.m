function [c, ceq] = jntrxncon_linForce(x, directionVector, LdirectionVector, maxAngle,LmaxAngle, A_f, LA_f, F_r0, LF_r0, GH_weight)
% This function, to be used within the optimizer FMINCON as the nonlinear
% constraint, returns whether the joint reaction constraint force lies
% inside a cone defined by maxAngle around the direction given by the
% vector generalDirection. It enforces effectively a directional constraint
% on the joint reaction force.
% Inputs are:
% * x: vector containing the activations of the actuators (could be muscle activations, controls, but also forces directly, if coherent throughout the whole process)
% * directionVector: vector giving the direction of optimal orientation for
%                    the joint reaction
% * maxAngle: maximum angle of accepted deviation wrt the optimal
%             orientation for the resulting joint reaction (defines a cone)
% * A_f: matrix containing the contribution of each actuator to the
%        reaction force
% * F_r0: value of the reaction force when x_i = 0, for all i
%* GH_weigth: weighting coefficient for the level of constraint violation.
%             If not passed, it is fixed to be 1
%
% The joint reaction force is found as a linear function of the muscle
% activations of the model considered, and takes the following formulation:
% force_vec = A_f * x' + F_r0;
% where the matrix A_f contains the difference between teh force caused by
% the activation of each actuator and the value F_r0 (found when all the
% actuators do not produce any force) and x is the vector of activations.
%
% The formulation implemented here assumes that the model state is defined
% and fixed for the given time interval. This justifies the linearized \
% formulation of the force
%
% author: Italo Belli (i.belli@tudelft.nl) 2022
%Edited By: Ibrahim Mohammed Hasan (imihasan@kth.se) 2023, KTH MoveAbility Lab,
%KTH Royal Institute of Technology, Stockholm, Sweden.

import org.opensim.modeling.*;

% set default weight for constraint violation if not passed by user
if nargin==5
    GH_weight = 1;
end

% computing the reaction force vector at the given joint
%Right
force_vec = A_f * x' + F_r0;
%Left
Lforce_vec = LA_f * x' + LF_r0;

% evaluate the relative angle between the reaction force and Vec_H2GC
%Right
cosTheta = max(min(dot(directionVector,force_vec)/(norm(directionVector)*norm(force_vec)),1),-1);
rel_angle = real(acosd(cosTheta));
%Left
LcosTheta = max(min(dot(LdirectionVector,Lforce_vec)/(norm(LdirectionVector)*norm(Lforce_vec)),1),-1);
Lrel_angle = real(acosd(LcosTheta));

% value of the constraint violation
c = [GH_weight*((rel_angle/maxAngle)^2 - 1); GH_weight*((Lrel_angle/LmaxAngle)^2 - 1)];             % direction must lie in a cone
ceq = [0; 0]; 
