% UAV Model with middle wing independently actuated i.e. left and right
clear all;
T = 10; % simulation time
dt= 0.002; % step time

N = T/dt; % number of data points in simulation
t=linspace(0,T,N);

% Gravity
g = 9.81;
% Air Density
rho = 1.225;

% Vehicle Properties

M_empty = 440;
M_max = 640;

N_jet = 36;
N_jet_front = 12;
N_jet_rear = 24;

v_max = 300; % [km/h]
v_cruise = 280; % [km/h]

% Assumptions
M_body = 400;
M_wings = M_max - M_body;
Thrust_max = M_max*g;
Ff_max = (Thrust_max*(N_jet_front/N_jet)/2); % Max magnitude for the input at the front wing
Fr_max = (Thrust_max*(N_jet_rear/N_jet)/2); % Max magnitude for the input at the back wing


% Relative Geometry
L = 2;

% Body properties -- Assume ellipsoid

x_b = L;
y_b = L/2;
z_b = L/2;

A_x = pi*y_b*z_b;
A_y = pi*x_b*z_b;


M = M_body; % Mass of fuselage

CDe = 0.10;
F0_max = 1/2*CDe*rho*(pi*y_b^2)*(v_max/3.6)^2;



% Rear wing dimensions -- Assume flat cuboid

x_wr = L/2;
y_wr = 4*L;
z_wr = L/10;
A_r = x_wr*y_wr;

% Front wing dimensions -- Assume flat cuboid

x_wf = L/2;
y_wf = 2*L;
z_wf = L/10;
A_f = x_wf*y_wf;

L_wing = 2*y_wf + 2*y_wr;

m = M_wings/ L_wing; % mass / unit length of wing [kg/m]

mf = y_wf*m;  % Mass of single front wing
mr = y_wr*m;  % Mass of single reA_r wing

% Moments of Inertia

% Body
Ixx = (1/5)*M*(y_b^2+z_b^2);
Iyy = (1/5)*M*(x_b^2+z_b^2);
Izz = (1/5)*M*(x_b^2+y_b^2);
I_body = diag([Ixx,Iyy,Izz]);

% Front Wing
Ixxwf = (1/12)*mf*(y_wf^2+z_wf^2);
Iyywf = (1/12)*mf*(x_wf^2+z_wf^2);
Izzwf = (1/12)*mf*(y_wf^2+x_wf^2);
I_front = diag([Ixxwf,Iyywf,Izzwf]);

% Rear Wing
Ixxwr = (1/12)*mr*(y_wr^2+z_wr^2);
Iyywr = (1/12)*mr*(x_wr^2+z_wr^2);
Izzwr = (1/12)*mr*(y_wr^2+x_wr^2);
I_rear = diag([Ixxwr,Iyywr,Izzwr]);

% Centre of mass locations
m_tot = M + 2*mr + 2*mf;
CoM = (M*x_b + 2*mf*2*x_b)/(m_tot);

% Geometry to wing centre of mass (same variable names as Maple)
Lf = y_b + y_wf/2; % Magnitude of centre of front wing in positive y-direction
Lr = y_b + y_wr/2; % Magnitude of centre of reA_r wing in positive y-direction
a = 2*x_b - CoM; % Magnitude of centre of front wing in positive x-direction
b = CoM; % Magnitude of centre of front wing in positive x-direction

% Wing damping and spring stiffness -- to return to neutral AoA.
cf = 0*mf;
kf = 0;
cr = 0*mr;
kr = 0;

Nd = 9; % degrees of freedom

% Flight mode
mode = 2; % 1 = vertical, 2 = horizontal
if mode == 1
    alpha0 = pi/2;
    alpha1 = pi/2;
    v0 = 0;
    ff0 = 1*b*m_tot*g/(2*(a+b));
    fr0 = 1*a*m_tot*g/(2*(a+b));
    
elseif mode == 2
    alpha1 = pi/10;
    alpha0 = 0;
    v0 = 50;
    
    %     [v0,alpha0] = SSFlight(v0,A_r,A_f,m_tot,rho,alpha1);
    %
    [Clf,Cdf] = AeroCoeffs(alpha0);
    [Clr,Cdr] = AeroCoeffs(alpha1);
    [Cll,Cdl] = AeroCoeffs(alpha1);
    
    flf = (1/2)*rho*v0^2*Clf*A_f; % front lift
    flr = (1/2)*rho*v0^2*Clr*A_r; %
    fll = (1/2)*rho*v0^2*Cll*A_r; %
    
    fdf = (1/2)*rho*v0^2*Cdf*A_f; % front drag
    fdr = (1/2)*rho*v0^2*Cdr*A_r; %
    fdl = (1/2)*rho*v0^2*Cdl*A_r; %
    fdx = (1/2)*rho*v0^2*CDe*A_x; %
    
    F0 = 2*fdf + 2*fdr + fdx;
    
    ff0 = 1/6*F0;
    fr0 = 1/3*F0;
    
end

% Note a negative z value indicates a positive altitude in the NED
% coordinate system.

q0 = [0;0;-1000;0;0;0;alpha0;alpha1;alpha1;zeros(Nd,1)];
q0(10) = v0;

q = q0;

Q(:,1) = q;

v_wind = 0;
theta_wind = 0;

% Disturbances

% Turn Wind noise on or off
Wind_noise = 0;

z0 = 4000;

for i = 2:N
    
    %     q(7:9) = q0(7:9);
    %     q(16:18) = [0;0;0];
    
    %     if q(3) <= .1
    %         q(3) = 0;
    %         if q(12) <= 0.1
    %             q(12) = 0;
    %         end
    %     end
    %
    
    %     theta_wind = theta_wind + Wind_noise*randn()/200;
    %     V = q(10:12);
    %     V_wind = (v_wind+Wind_noise*randn()/1000)*[cos(theta_wind);sin(theta_wind);0];
    %     v_rel = V+V_wind;
    %
    %     alpha_rel = atan2(v_rel(3),norm(v_rel(1:2)));
    %     psi_head = atan2(v_rel(2),v_rel(1));
    %     psi_rel = psi_head - q(6);
    %
    %     v_bx = norm(v_rel)*cos(psi_rel);
    %     v_by = norm(v_rel)*sin(psi_rel);
    
    theta_wind = theta_wind + Wind_noise*randn()/200;
    V = q(10);
    v_rel = V;
    
    alpha_rel = atan2(0,norm(v_rel(1)));
    psi_head = atan2(0,v_rel(1));
    psi_rel = psi_head - q(6);
    
    v_bx = norm(v_rel)*cos(psi_rel);
    v_by = norm(v_rel)*sin(psi_rel);
    
    
    %     figure(2)
    %     plot(theta_wind,psi_rel,'r.',theta_wind,v_bx,'b.')
    %     hold on
    G = RiemannianMetric(q,Lf,Lr,a,b,M,mf,mr,I_body,I_front,I_rear);
    
    % Control Forces
    fc = - 10* (q(10)-v0);
    
    ffl = ff0 + fc;
    ffr = ff0 + fc;
    frl = fr0 + fc;
    frr = fr0 + fc;
    
    T_on = 1;
    k_t  = 0.01;
    Tf   = T_on*(k_t*((q(3)-q0(3))-.0*(q(16)-0)));
    Trl  = T_on*(k_t*((q(3)-q0(3))-.0*(q(17)-0)));
    Trr  = T_on*(k_t*((q(3)-q0(3))-.0*(q(18)-0)));
    
    u = [ffl;ffr;frl;frr;Tf;Trl;Trr];
    
    % Aerodynamic Forces
    
    [Clf,Cdf] = AeroCoeffs(q(7)-alpha_rel);
    [Cll,Cdl] = AeroCoeffs(q(8)-alpha_rel);
    [Clr,Cdr] = AeroCoeffs(q(9)-alpha_rel);
    
    % Lift forces
    
    fl_on = 1;
    
    flf = fl_on*(1/2)*rho*v_bx^2*Clf*A_f; % front lift
    fll = fl_on*(1/2)*rho*v_bx^2*Cll*A_r; % reA_r left lift
    flr = fl_on*(1/2)*rho*v_bx^2*Clr*A_r; % reA_r right lift
    
    % Drag forces
    
    fd_on = 1*sign(v_bx);
    
    fdf = fd_on*(1/2)*rho*v_bx^2*Cdf*A_f; % front lift
    fdl = fd_on*(1/2)*rho*v_bx^2*Cdl*A_r; % reA_r left lift
    fdr = fd_on*(1/2)*rho*v_bx^2*Cdr*A_r; % reA_r right lift
    
    fbx = fd_on*(1/2)*rho*v_bx^2*CDe*A_x;
    fby = 1*(1/2)*rho*v_by^2*.013*A_y;
    
    f_aero = [flf;fll;flr;fdf;fdl;fdr;fbx;fby];
    
    % Calculate contribution of each force
    F_drift = DriftForces(q,Lf,Lr,a,b,m_tot,mf,mr,I_body,I_front,I_rear);
    F_aero = AerodynamicForces(q,f_aero,Lf,Lr,a,b);
    F_actuated = ActuatedForces(q,u,Lf,Lr,a,b);
    
    qdot = [q((Nd+1):(2*Nd));...
        inv(G)*(F_drift + F_aero + F_actuated - [zeros(1,Nd-3),cf,cr,cr]*[zeros(Nd-3,1);q((2*Nd-2):2*Nd)] - [zeros(1,Nd-3),kf,kr,kr]*[zeros(Nd-3,1);q((Nd-2):Nd)])];
    
    q = q + qdot*dt;
    
    Q(:,i) = q;
    F(:,i) = f_aero;
    U(:,i) = u;
    ALPHA(:,i) = alpha_rel;
    PSI(:,i) = psi_rel;
    
    V_REL(:,i) = [v_bx;v_by];
    
end

figure(1)
set(gcf,'color','w');
clf
subplot(3,2,[1])
plot(t,Q(1:2,:),t,-Q(3,:),t,Q(10:11,:),t,-Q(12,:))
grid on
legend({'$x$','$y$','$z$','$\dot{x}$','$\dot{y}$','$\dot{z}$'},'Interpreter','latex','Location','best')
subplot(3,2,[3])
plot(t,Q(4:6,:),t,Q(13:15,:))
grid on

legend({'$\theta$','$\phi$','$\psi$','$\dot{\theta}$','$\dot{\phi}$','$\dot{\psi}$'},'Interpreter','latex','Location','best')
subplot(3,2,5)
plot(t,Q(7:9,:),t,Q(16:18,:))
grid on
legend({'$\alpha$','$\alpha_l$','$\alpha_r$','$\dot{\alpha}$','$\dot{\alpha_l}$','$\dot{\alpha_r}$'},'Interpreter','latex','Location','best')
subplot(3,2,[4])
plot(t,F)
grid on
legend({'$F_{l_f}$','$F_{l_l}$','$F_{l_r}$','$F_{d_f}$','$F_{d_l}$','$F_{d_r}$'},'Interpreter','latex','Location','best')
subplot(3,2,6)
plot(t,U)
grid on
legend({'$F_{f_l}$','$F_{f_r}$','$F_{r_l}$','$F_{r_r}$','$T_{f}$','$T_{r_l}$','$T_{r_r}$'},'Interpreter','latex','Location','best')
subplot(3,2,2)
plot(t,ALPHA,t,PSI,t,V_REL)
grid on
legend({'$\alpha_{rel}$','$\psi_{rel}$','$v_{b_x}$','$v_{b_y}$'},'Interpreter','latex','Location','best')

