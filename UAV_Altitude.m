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
k = 0.99*m0;

Nd = 2; % degrees of freedom

rho = 1.225;

T1 = 0*(sin(.5*t)+0);
T2 = 0*(sin(.5*t)+0);

alpha0=pi/20;
Cl0 = 3*(sin(alpha0+pi/40));
Cd0 = 10*(sin(alpha0-pi/2)+1);
v0 = sqrt(m0*g*cos(alpha0)/(Cl0*cos(alpha0)-Cd0*sin(alpha0)));
F0 = (1/2)*(Cl0*sin(alpha0)+Cd0*cos(alpha0))*v0^2-m0*g*sin(alpha0)/2;

q0 = [0;0;0;0];

q = q0;

Q(:,1) = q;

theta0 = 0;

z0 = 1000*t./t;

for i = 2:N
    
    v = v0 + 1*randn()/1;
    
    if q(2) > pi/6
        Cll = 3*(sin(pi/6+pi/40));
        Cdl = 1*10*(sin(pi/6-pi/2)+1);
        q(2) = pi/6;
        if q(4) > 0
            q(4) = 0;
        end
    elseif q(2) < -pi/6
        Cll = 3*(sin(-pi/6+pi/40));
        Cdl = 1*10*(sin(-pi/6-pi/2)+1);
        q(2) = -pi/6;
        if q(4) < 0
            q(4) = 0;
        end
    else
        Cll = 3*(sin(q(2)+pi/40));
        Cdl = 1*10*(sin(q(2)-pi/2)+1);
    end
    
    G = [
        M + 2*m, 0;
        0, 2*Iyyw];
    
    % Actual Forces
    fll = 1*(1/2)*rho*v^2*Cll;
    flr = fll;
    fdl = 1*sign(v)*(1/2)*rho*v^2*Cdl;
    fdr = fdl;
    
    T1 = -0.5*[0.2 0.9]*[(q(1)-z0(i));(q(3)-0)];
    
    GammaPrime = ...
        [
        -(M+2*m)*g + fll + flr;
        T1-c*q(4)-k*(q(2)-alpha0)];...
        
    qdot = [q((Nd+1):(2*Nd));...
        inv(G)*(GammaPrime)];
    
    q = q + qdot*dt;
    
    
    Q(:,i) = q;
    F(:,i) = [fll;fdl];
    CD(:,i) = [Cll;Cdl];
    U(:,i) = [T1];
    V(:,i) = [v];
    
end

figure(1)
clf
set(gcf,'color','w');

subplot(3,2,1)
plot(t,Q(1,:))
grid on
legend({'$z$'},'Interpreter','latex','Location','best')
subplot(3,2,3)
plot(t,Q(3,:))
grid on

legend({'$\dot{z}$'},'Interpreter','latex','Location','best')
subplot(3,2,5)
plot(t,Q(2,:))
grid on

legend({'$\alpha$'},'Interpreter','latex','Location','best')
subplot(3,2,[2])
plot(t,F)
grid on

legend({'$F_{l}$','$F_{d}$'},'Interpreter','latex','Location','best')
subplot(3,2,4)
plot(t,CD)
grid on

legend({'$C_{l}$','$C_{d}$'},'Interpreter','latex','Location','best')
subplot(3,2,6)
plot(t,U)
grid on

legend({'$T$'},'Interpreter','latex','Location','best')

figure(2)
set(gcf,'color','w');

clf
plot(t,V)
grid on

