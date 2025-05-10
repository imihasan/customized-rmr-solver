function c = cost(x, w, Vec_H2GC, Vec_GC2GE, A_f, F_r0, tag_cost)

% Author: Ibrahim Mohammed I. Hasan (imihasan@kth.se) 2025

% computing the reaction force vector at the given joint
force_vec = A_f * x' + F_r0;


% evaluate the relative angle between the reaction force and Vec_H2GC
cosTheta = max(min(dot(Vec_H2GC,force_vec)/(norm(Vec_H2GC)*norm(force_vec)),1),-1);
%cosTheta = dot(Vec_H2GC,force_vec)/(norm(Vec_H2GC)*norm(force_vec));
theta = real(acosd(cosTheta));

%Obtain the rotation matrix
norm_Vec_H2GC=Vec_H2GC/norm(Vec_H2GC);
xhat=-norm_Vec_H2GC;
norm_Vec_GC2GEE=Vec_GC2GE/norm(Vec_GC2GE);
Vxdot=dot(norm_Vec_GC2GEE,xhat);
Vy=norm_Vec_GC2GEE-Vxdot*xhat;
yhat=Vy/norm(Vy);
zhat=cross(xhat,yhat);
R=[xhat; yhat; zhat];


fv_rotated = R*force_vec;

if tag_cost=="conditional"
    ratio=(fv_rotated(2)^2+fv_rotated(3)^2)/fv_rotated(1)^2;
    indx=ratio;
    if indx>=(0.5)^2
        Cond_Penalty=((indx));
    else
        Cond_Penalty=0;
    end
    c = sum(w.*(x.^2))+3*Cond_Penalty; %Conditional Penalty
elseif tag_cost=="planar"
    ratio=(fv_rotated(2)^2+fv_rotated(3)^2)/fv_rotated(1)^2;
    planar_Penalty=ratio;
    c = sum(w.*(x.^2))+3*planar_Penalty; %Planar Penalty
elseif tag_cost=="curve"
    rc=0.1; %raduis of the glenoud cavity sphere
    dch=0.09-norm(Vec_H2GC);
    phl=(2*dch*cosd(180-theta)+sqrt((2*dch*cosd(180-theta))^2-4*(dch^2-rc^2)))/2;
    pgl_prime=phl*abs(sind(theta));
    beta=asind(pgl_prime/rc);
    alpha=theta-beta;
    fs=norm(fv_rotated)*sind(alpha); %tangential component to the curvature
    fc=norm(fv_rotated)*cosd(alpha); %nomal component to the curvature
    curve_Penalty=fs^2/fc^2;
    c = sum(w.*(x.^2))+10*curve_Penalty; %Curve Penalty
elseif tag_cost=="no_penalty"
    c = sum(w.*(x.^2)); %no penalty
end

end