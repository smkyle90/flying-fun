% ACTUATED AIRFOIL
clear all;
T = 50; % simulation time
dt= 0.001; % step time

N = T/dt; % number of data points in simulation
t=linspace(0,T,N);

% Gravity
g = 9.81;
% Air Density
rho = 1.225;

% Relative Geometry
L = 0.5;
W = 8*L;
m = 100;
k=0;
c=0;
a = -0.75;
I = (1/12)*m*(2*L)^2;

Nd = 3; % degrees of freedom

v_wind = 0;
theta_max = pi/6;

y0 = 3000;
theta0 = pi/20;

[Cl0,Cd0] = AeroCoeffs(theta0);

f0 = (Cd0/Cl0)*m*g;
v0 = sqrt(2*m*g/(rho*Cl0*2*L*W))+v_wind;

fl0 = (1/2)*rho*(v0)^2*Cl0*(2*L*W); % front lift
fd0 = (1/2)*rho*(v0)^2*Cd0*(2*L*W); % front lift

f0 = fd0;
t0 = L*sin(theta0)*f0;

q0 = [0;2000;theta0;0;0;0];

q = q0;
Q(:,1) = q;

for i = 2:N
    
    v_wind = v_wind+0*randn()/50;
    
    %     figure(2)
    %     plot(theta_wind,psi_rel,'r.',theta_wind,v_bx,'b.')
    %     hold on
    G = [m, 0, m*L*sin(q(3));...
        0, m,-m*L*cos(q(3));...
        m*L*sin(q(3)), -m*L*cos(q(3)), I + m*L^2];
    % Aerodynamic Forces
    
    [Cl,Cd] = AeroCoeffs(q(3));
    
    v_air = q(4)-v_wind;
    % Lift forces
    fl_on = 1;
    fl = fl_on*(1/2)*rho*(v_air)^2*Cl*(2*L*W); % front lift
    
    % Drag forces
    fd_on = 1*sign(v_air);
    fd = fd_on*(1/2)*rho*(v_air)^2*Cd*(2*L*W); % front lift
    
    X_drift = ...
        [-m*L*q(6)^2*cos(q(3));...
        -m*L*q(6)^2*sin(q(3))-m*g;...
        m*g*L*cos(q(3))];
    
    dV = [q(1)-v0*t(i);0;(q(3)-(-(theta_max)*atan((q(2)-y0)/200)/(pi/pi)))];
    
    Fd = a*q(4:6);
    
    f1 = min(max(0,(m)*(-1*dV(1) + Fd(1)) + f0),1000);
    t1 = (4*m)*(-1*dV(3) + 1*Fd(3) - .01*q(5)) + t0;
    
    X_aero = ...
        [0, -1;...
        1,  0;...
        -L*cos(q(3)),-L*sin(q(3))];
    
    X_force = ...
        [1, 0;...
        0, 0;...
        0, 1];
    
    f_aero = [fl;fd];
    u = [f1;t1];
    
    qdot = [q(Nd+1:2*Nd);...
        inv(G)*(X_drift + X_aero*f_aero + X_force*u)];
    
    
    q = q + qdot*dt;
    
    if q(2) < 0.1
        q(2) = 0;
        if q(5) <= 0
            q(5)=0;
        end
    end
    
    Q(:,i) = q;
    F(:,i) = [fl;fd];
    %     CD(:,i) = [Cll;Clr;Cdl;Cdr];
    U(:,i) = u;
    V(:,i) = [v_wind];
    
end

figure(1)
set(gcf,'color','w');
clf
subplot(3,2,[1])
plot(t,Q(3,:),t,Q(6,:))
hold on
plot(t,theta0*ones(1,N),'k--',t,theta_max*ones(1,N),'b--',t,-theta_max*ones(1,N),'b--')
grid on
legend({'$\theta$','$\dot{\theta}$'},'Interpreter','latex','Location','best')
subplot(3,2,[3])
plot(t,Q(1,:),t,Q(4,:))
grid on
legend({'$x$','$\dot{x}$'},'Interpreter','latex','Location','best')
subplot(3,2,[5])
plot(t,Q(2,:),t,Q(5,:))
grid on
legend({'$y$','$\dot{y}$'},'Interpreter','latex','Location','best')
subplot(3,2,[2])
plot(t,F)
grid on
legend({'$F_{l}$','$F_{d}$'},'Interpreter','latex','Location','best')
subplot(3,2,[4])
plot(t,U)
grid on
legend({'$F$','$T$'},'Interpreter','latex','Location','best')
subplot(3,2,[6])
plot(t,V,t,Q(4,:),t,Q(4,:)-V)
grid on
legend({'$V_{w}$','$V_{g}$','$V_{a}$'},'Interpreter','latex','Location','best')

