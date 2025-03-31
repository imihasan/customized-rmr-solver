function [c, ceq] = jntrxncon_linForce(x, Vec_H2GC, Vec_GC2GE, A_f, F_r0,tag_con)
% This function, to be used within the optimizer FMINCON as the nonlinear
% constraint, returns whether the joint reaction constraint force lies
% inside a cone defined by maxAngle around the direction given by the
% vector generalDirection. It enforces effectively a directional constraint
% on the joint reaction force.
% Inputs are:
% * x: vector containing the activations of the actuators (could be muscle activations, controls, but also forces directly, if coherent throughout the whole process)
%
% * directionVector: vector giving the direction of optimal orientation for
%                    the joint reaction
%
% * A_f: matrix containing the contribution of each actuator to the
%
%        reaction force
% * F_r0: value of the reaction force when x_i = 0, for all i
%
% * tag_con: string to define which constraint to be enforced: "point",
%            "circle", "ellipse" or "polynomial"
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

import org.opensim.modeling.*;

% set default weight for constraint violation if not passed by user
% if nargin==5
%     GH_weight = 1;
% end

% computing the reaction force vector at the given joint
force_vec = A_f * x' + F_r0;

% evaluate the relative angle between the reaction force and Vec_H2GC
cosTheta = max(min(dot(Vec_H2GC,force_vec)/(norm(Vec_H2GC)*norm(force_vec)),1),-1);
%cosTheta = dot(Vec_H2GC,force_vec)/(norm(Vec_H2GC)*norm(force_vec));
rel_angle = real(acosd(cosTheta));


%These compuations of the rotation matrix were prone to gimbal lock
%another derivation for the rotation matrix is introduced below
% beta_angle = atan(norm_Vec_H2GC(3)/norm_Vec_H2GC(1));
% alpha_angle = atan(norm_Vec_H2GC(3)/(sin(beta_angle)*norm_Vec_H2GC(2)));
% Ry = [cos(beta_angle) 0 sin(beta_angle); 0 1 0; -sin(beta_angle) 0 cos(beta_angle)];
% Rz = [cos(alpha_angle) -sin(alpha_angle) 0; sin(alpha_angle) cos(alpha_angle) 0; 0 0 1];

%A new derivation of the roation matrix that transforms to th glenoid
%cavity frame of reference
norm_Vec_H2GC=Vec_H2GC/norm(Vec_H2GC);
xhat=-norm_Vec_H2GC;
norm_Vec_GC2GEE=Vec_GC2GE/norm(Vec_GC2GE);
Vxdot=dot(norm_Vec_GC2GEE,xhat);
Vy=norm_Vec_GC2GEE-Vxdot*xhat;
yhat=Vy/norm(Vy);
zhat=cross(xhat,yhat);
R=[xhat; yhat; zhat]; %The rotation matrix that represents vectors in the glenoid frame

fv=R*force_vec;
norm_fv_in_ground = force_vec/norm(force_vec);
norm_fv_rotated = R*norm_fv_in_ground*(norm(1000*Vec_H2GC)/cosTheta);




%%


% value of the constraint violation

if tag_con=="polynomial"
    %Polynomial formulation
    if -norm_fv_rotated(2)>=0
        theta=atan2d(-norm_fv_rotated(2),-norm_fv_rotated(3));
    else
        theta=atan2d(-norm_fv_rotated(2),-norm_fv_rotated(3))+360;
    end
    L=(-5.430378773797311e-12*theta^6+5.880441464447251e-09*theta^5+-2.309728860559176e-06*theta^4+3.899916118022262e-04*theta^3+-0.024526120335986*theta^2+0.116047310713347*theta+50.915928515941545)/100;
    c=(fv(2)^2+fv(3)^2)/fv(1)^2-L^2; %polynomial Constraint

elseif tag_con=="circle"
    h2gc=R*Vec_H2GC'*1000; %vector from the glenohumeral head to the gleoid
    % center convertedto mm
    radius = norm(h2gc)*0.5; %This establishes a circular border corresponds to
    % 0.5 stability ratio. 

    %This was the original formulation in the rmr solver to enforce the circle
    %constraint
    %c = ((rel_angle/maxAngle)^2 - 1);% direction must lie in a cone
    %The above constraint was replaced by this one
    c=(norm([-norm_fv_rotated(2) -norm_fv_rotated(3)])/radius)^2-1; %circle constraint

elseif tag_con=="point"
    c=rel_angle; %point constraint
elseif tag_con=="ellipse"
    c=[(fv(3)^2/0.61^2+fv(2)^2/0.34^2)/fv(1)^2-1]; %ellipse constraint
end

ceq = 0; 