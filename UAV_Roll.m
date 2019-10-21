% UAV Model with middle wing independently actuated i.e. left and right
clear all;
T = 100; % simulation time
dt= 0.01; % step time

N = T/dt; % number of data points in simulation
t=linspace(0,T,N);

M = 10;
m = 5;
m0 = M + 2*m;

L = 2;

g = 9.81;

% Body properties -- Assume ellipsoid

x_b = 2*L;
y_b = L;
z_b = L;

Ixx = (1/5)*M*(y_b^2+z_b^2);
Iyy = (1/5)*M*(x_b^2+z_b^2);
Izz = (1/5)*M*(x_b^2+y_b^2);

% Wing properties -- Assume flat cuboid

x_w = L/2;
y_w = 2*L;
z_w = L/10;

Ixxw = (1/12)*m*(y_w^2+z_w^2);
Iyyw = (1/12)*m*(x_w^2+z_w^2);
Izzw = (1/12)*m*(x_w^2+y_w^2);

c = 0.5*m0;
k = 1*m0;

Nd = 3; % degrees of freedom

rho = 1.225;

T1 = 0*(sin(.5*t)+0);
T2 = 0*(sin(.5*t)+0);

alpha0=0;
Cl0 = 3*(sin(alpha0+pi/40));
Cd0 = 10*(sin(alpha0-pi/2)+1);
v0 = sqrt(m0*g*cos(alpha0)/(Cl0*cos(alpha0)-Cd0*sin(alpha0)));
F0 = (1/2)*(Cl0*sin(alpha0)+Cd0*cos(alpha0))*v0^2-m0*g*sin(alpha0)/2;

q0 = [0;pi/100;-pi/40;0.0000000000;0.0000000;0.0000000];

q = q0;

Q(:,1) = q;

theta0 = 0*t;
omega0 = 0*t;

for i = 2:N
    
    v = v0;
    
    vl = v + 1*randn()/3;
    vr = v + 1*randn()/3;
    
%     if q(1) > pi/4
%         q(1) = pi/4;
%         if q(4) > 0
%             q(4) = 0;
%         end
%         
%     elseif q(1) < -pi/4
%         q(1) = -pi/4;
%         if q(4) < 0
%             q(4) = 0;
%         end
%     end
    
    if q(2) > pi/6
        Cll = 3*(sin(pi/6+pi/40));
        Cdl = 1*10*(sin(pi/6-pi/2)+1);
        q(2) = pi/6;
        if q(5) > 0
            q(5) = 0;
        end
    elseif q(2) < -pi/6
        Cll = 3*(sin(-pi/6+pi/40));
        Cdl = 1*10*(sin(-pi/6-pi/2)+1);
        q(2) = -pi/6;
        if q(5) < 0
            q(5) = 0;
        end
    else
        Cll = 3*(sin(q(2)+pi/40));
        Cdl = 1*10*(sin(q(2)-pi/2)+1);
        
    end
    
    if q(3) > pi/6
        Clr = 3*(sin(pi/6+pi/40));
        Cdr = 1*10*(sin(pi/6-pi/2)+1);
        q(3) = pi/6;
        if q(6) > 0
            q(6) = 0;
        end
    elseif q(3) < -pi/6
        Clr = 3*(sin(-pi/6+pi/40));
        Cdr = 1*10*(sin(-pi/6-pi/2)+1);
        q(3) = -pi/6;
        if q(6) < 0
            q(6) = 0;
        end
    else
        Clr = 3*(sin(q(3)+pi/40));
        Cdr = 1*10*(sin(q(3)-pi/2)+1);
    end
    
    % Modelled Forces
    fll = 1*(1/2)*rho*v^2*Cll;
    flr = 1*(1/2)*rho*v^2*Clr;
    
    fdl = 1*sign(vl)*(1/2)*rho*v^2*Cdl;
    fdr = 1*sign(vr)*(1/2)*rho*v^2*Cdr;
        
    T1 = -[.1 .5]*[(q(1)-theta0(i));(q(4)-omega0(i))];
    
    Tl = T1;
    Tr = -T1;
%     T1 = L*sin(q(1))*(sin(q(1))*(L^2*m - Iyyw + Izzw)*(q(5)^2 + q(6)^2)*((sin(q(2))*fll - sin(q(3))*flr)*cos(q(1)) + cos(q(2))*fdl - cos(q(3))*fdr)*cos(q(1)) + L*((sin(q(2))*fll - sin(q(3))*flr)*(cos(q(2))*fll - cos(q(3))*flr - fll + flr)*cos(q(1))^3 + (2*cos(q(2))^2*fdl*fll + ((-fdl*flr - fdr*fll)*cos(q(3)) - fdl*(fll - flr))*cos(q(2)) + 2*cos(q(3))^2*fdr*flr + fdr*(fll - flr)*cos(q(3)) + sin(q(3))*(fdl*flr + fdr*fll)*sin(q(2)) - fdl*fll - fdr*flr)*cos(q(1))^2 + ((-sin(q(2))*fdl^2 + sin(q(3))*fdl*fdr)*cos(q(2)) + fdr*(sin(q(2))*fdl - sin(q(3))*fdr)*cos(q(3)) + (fll - flr)*(sin(q(2))*fll - sin(q(3))*flr))*cos(q(1)) + (fll - flr)*(cos(q(2))*fdl - cos(q(3))*fdr)))*k_c*(q(1) - theta0 + q(4) - 0)/(2*(q(5)^2 + q(6)^2)^2*(L^2*m - Iyyw + Izzw)^2*cos(q(1))^4 - 2*(q(5)^2 + q(6)^2)^2*(L^2*m - Iyyw + Izzw)^2*cos(q(1))^2 - 4*sin(q(1))*(L^2*m - Iyyw + Izzw)*L*(q(5)^2 + q(6)^2)*((cos(q(2))*fll - cos(q(3))*flr - fll + flr)*cos(q(1))^2 + (-sin(q(2))*fdl + sin(q(3))*fdr)*cos(q(1)) + fll - flr)*cos(q(1)) - 3*((cos(q(2))^2*fll^2 - (4*fll*(cos(q(3))*flr + fll - flr)*cos(q(2)))/3 + cos(q(3))^2*flr^2 + (4*flr*(fll - flr)*cos(q(3)))/3 + (2*sin(q(2))*sin(q(3))*fll*flr)/3 - (4*fll*flr)/3 + fll^2/3 + flr^2/3)*cos(q(1))^4 + ((-2*sin(q(2))*fdl*fll + (2*sin(q(3))*(fdl*flr + 2*fdr*fll))/3)*cos(q(2)) + (((4*fdl*flr)/3 + (2*fdr*fll)/3)*sin(q(2)) - 2*sin(q(3))*fdr*flr)*cos(q(3)) + 4*(fll - flr)*(sin(q(2))*fdl - sin(q(3))*fdr)/3)*cos(q(1))^3 + ((-fdl^2 - fll^2/3)*cos(q(2))^2 + ((2*cos(q(3))*fdl*fdr)/3 + (4*fll*(fll - flr))/3)*cos(q(2)) + (-fdr^2 - flr^2/3)*cos(q(3))^2 - (4*flr*(fll - flr)*cos(q(3)))/3 - 4*(fdl*fdr + fll*flr/2)*sin(q(3))*sin(q(2))/3 + (8*fll*flr)/3 + (2*fdl^2)/3 + (2*fdr^2)/3 - fll^2 - flr^2)*cos(q(1))^2 + ((2*fdl*(sin(q(2))*fll - sin(q(3))*flr)*cos(q(2)))/3 - (2*fdr*(sin(q(2))*fll - sin(q(3))*flr)*cos(q(3)))/3 - 4*(fll - flr)*(sin(q(2))*fdl - sin(q(3))*fdr)/3)*cos(q(1)) + cos(q(2))^2*fdl^2/3 - (2*cos(q(2))*cos(q(3))*fdl*fdr)/3 + cos(q(3))^2*fdr^2/3 + (2*(fll - flr)^2)/3)*L^2);
    
    G = [
        2*L^2*m + Ixx + 2*Ixxw, 0, 0;...
        0, (Iyyw - Izzw)*cos(q(1))^2 + Izzw, 0;...
        0, 0, (Iyyw - Izzw)*cos(q(1))^2 + Izzw];
    
    % Actual Forces
    fll = 1*(1/2)*rho*vl^2*Cll;
    flr = 1*(1/2)*rho*vr^2*Clr;
    
    fdl = 0*sign(vl)*(1/2)*rho*vl^2*Cdl;
    fdr = 0*sign(vr)*(1/2)*rho*vr^2*Cdr;
    
    GammaPrime = ...
        [-cos(q(1))*sin(q(1))*(q(5)^2 + q(6)^2)*(Iyyw - Izzw) + L*(fll - flr);
         Tl;...
         Tr];...
        
    qdot = [q((Nd+1):(2*Nd));...
        inv(G)*(GammaPrime - [zeros(Nd-2,1);c*q(2*Nd-1);c*q(2*Nd)] - [zeros(Nd-2,1);k*q(Nd-1);k*q(Nd)])];
    
    q = q + qdot*dt;
    
    
    Q(:,i) = q;
    F(:,i) = [fll;flr;fdl;fdr];
    CD(:,i) = [Cll;Clr;Cdl;Cdr];
    U(:,i) = [Tl;Tr];
    V(:,i) = [vl;vr];
    
end

figure(1)
clf
subplot(3,2,1)
plot(t,Q(1,:))
grid on
legend({'$\theta$'},'Interpreter','latex','Location','best')
subplot(3,2,3)
plot(t,Q(4,:))
legend({'$\dot{\theta}$'},'Interpreter','latex','Location','best')
subplot(3,2,5)
plot(t,Q(2:3,:))
% hold on
% plot(Q(5:6,:))
legend({'$\alpha_l$','$\alpha_r$','$\dot{\alpha_l}$','$\dot{\alpha_r}$'},'Interpreter','latex','Location','best')
subplot(3,2,[2])
plot(t,F)
legend({'$F_{l_l}$','$F_{l_r}$','$F_{d_l}$','$F_{d_r}$'},'Interpreter','latex','Location','best')
subplot(3,2,4)
plot(t,CD)
legend({'$C_{l_l}$','$C_{l_r}$','$C_{d_l}$','$C_{d_r}$'},'Interpreter','latex','Location','best')
subplot(3,2,6)
plot(t,U)
legend({'$T_{l}$','$T_{r}$'},'Interpreter','latex','Location','best')

figure(2)
clf
plot(t,V)
