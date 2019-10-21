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

q0 = [0;0;0;v0;0;0];

q = q0;

Q(:,1) = q;

theta0 = 0;

z0 = 1200+t;
x0 = v0*t;

for i = 2:N
    
    v = q(4);
    if q(2) < 0
        q(2) = 0;
        if q(5) < 0
            q(5) = 0;
        end
    end
    
    if q(3) > pi/6
        Cl = 3*(sin(pi/6+pi/40));
        Cd = 1*10*(sin(pi/6-pi/2)+1);
        q(3) = pi/6;
        if q(6) > 0
            q(6) = 0;
        end
    elseif q(3) < -pi/6
        Cl = 3*(sin(-pi/6+pi/40));
        Cd = 1*10*(sin(-pi/6-pi/2)+1);
        q(3) = -pi/6;
        if q(6) < 0
            q(6) = 0;
        end
    else
        Cl = 3*(sin(q(3)+pi/40));
        Cd = 1*10*(sin(q(3)-pi/2)+1);
    end
    
    G = [
        M + 2*m, 0,0;
        0, M+2*m,0;
        0,0, 2*Iyyw];
    
    % Actual Forces
    fl = 1*(1/2)*rho*v^2*Cl;
    fd = 1*sign(v)*(1/2)*rho*v^2*Cd;
    
    f = -2*(q(1)-x0(i))-2*(q(4)-v0);
    
    if f < 0
        f=0;
    end
    T1 = -2*[0.005,0.05]*[(q(2)-z0(i));(q(5)-0)];
    GammaPrime = ...
        [-2*fd + 2*cos(q(3))*f;
        -m0*g + 2*fl - 2*sin(q(3))*f;
        T1-c*q(6)-k*(q(3)-alpha0)];...
        
    qdot = [q((Nd+1):(2*Nd));...
        inv(G)*(GammaPrime)];
    
    q = q + qdot*dt;
    
    
    Q(:,i) = q;
    F(:,i) = [fl;fd];
    CD(:,i) = [Cl;Cd];
    U(:,i) = [f,T1];
    V(:,i) = [v];
    
end

figure(1)
set(gcf,'color','w');
clf
subplot(3,2,1)
plot(t,Q(1:2,:))
grid on
hold on
plot(t,z0,t,x0)

legend({'$x$','$z$','$z_0$','$x_0$'},'Interpreter','latex','Location','best')
subplot(3,2,3)
plot(t,Q(4:5,:))
grid on

legend({'$\dot{x}$','$\dot{z}$'},'Interpreter','latex','Location','best')
subplot(3,2,5)
plot(t,Q(3,:))
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

legend({'$f$','$T$'},'Interpreter','latex','Location','best')

figure(2)
set(gcf,'color','w');
clf
plot(t,V)
grid on

