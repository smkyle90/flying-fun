% This function determines the  effect of actuated forces on the
% body.

function F_actuated = ActuatedForces(q,U,Lf,Lr,a,b)

theta = q(4);
phi = q(5);
psi = q(6);
alpha = q(7);
alpha_l = q(8);
alpha_r = q(9);

F_actuated = ...
    [cos(psi)*cos(phi + alpha), cos(psi)*cos(phi + alpha), cos(psi)*cos(phi + alpha_l), cos(psi)*cos(phi + alpha_r),0,0,0;...
    sin(psi)*cos(phi + alpha), sin(psi)*cos(phi + alpha), sin(psi)*cos(phi + alpha_l), sin(psi)*cos(phi + alpha_r),0,0,0;...
    -sin(phi + alpha), -sin(phi + alpha), -sin(phi + alpha_l), -sin(phi + alpha_r),0,0,0;...
    Lf*sin(alpha)*cos(theta), -Lf*sin(alpha)*cos(theta), Lr*sin(alpha_l)*cos(theta), -Lr*sin(alpha_r)*cos(theta),0,0,0;...
    -Lf*cos(alpha)*sin(theta) + a*sin(alpha), Lf*cos(alpha)*sin(theta) + a*sin(alpha), -Lr*cos(alpha_l)*sin(theta) - b*sin(alpha_l), Lr*cos(alpha_r)*sin(theta) - sin(alpha_r)*b,0,0,0;...
    cos(theta)*Lf*cos(phi + alpha), -cos(theta)*Lf*cos(phi + alpha), Lr*cos(phi + alpha_l)*cos(theta), -Lr*cos(phi + alpha_r)*cos(theta),0,0,0;...
    0, 0, 0, 0, 1, 0, 0;...
    0, 0, 0, 0, 0, 1, 0;...
    0, 0, 0, 0, 0, 0, 1]*U;


