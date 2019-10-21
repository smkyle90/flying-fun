% This function determines the  coefficient of lift and drag for an
% airfoill at an angle of attack equal to alpha

function [Cl,Cd] = AeroCoeffs(alpha)

% Lift coefficient

if alpha > pi/6
    Cl = 0*(sin(pi/6+pi/40));
elseif alpha < -pi/6
    Cl = 0*(sin(-pi/6+pi/40));
else
    Cl = 3*sin(alpha+pi/40);
end

% Drag coefficient

Cd = min(2,1*5*(sin(alpha-pi/2)+1));
