%% Parameter Initializations
clear
close all
fixed_point = [0,1; 0.5,1.5;1,1]';
max(size(fixed_point));

m = 1; %kg
g = 9.81 ; % gravity 

%% Initial condition
p_0 = [0.5,0]';
tau_min = 1;
tau_max = 30;

%% Initial Configuration
figure(1)
plot(p_0(1),p_0(2),'ko'); 
hold on
plot(fixed_point(1,:),fixed_point(2,:), 'ro')
hold on
X_0 = [fixed_point(1,:); p_0(1)*ones(1,max(size(fixed_point)))];
Y_0 = [fixed_point(2,:); p_0(2)*ones(1,max(size(fixed_point)))];
plot(X_0,Y_0)
title('Scheme of the three cables system');
grid on ; xlabel('Time (s)') ;

%% Tension distribution of initial condition
U = fixed_point - p_0*ones(1, max(size(fixed_point)));
W_0 = normc(U);
G = [0,-m*g]';
N_0 = null(W_0);
x0s = W_0\-G;
lb = tau_min*ones(max(size(fixed_point)),1) - x0s;
ub = tau_max*ones(max(size(fixed_point)),1) - x0s;
fun = @(x) norm(x0s+ N_0*x);
cond_A = [N_0;-N_0];
cond_b = [ub; -lb];
lamda = fmincon(fun,-3.8,cond_A,cond_b);
Tau = x0s+ N_0*lamda;

fprintf('The Initial cable tension distribution is: [');
fprintf('%g, ', Tau(1:end-1));
fprintf('%g]\n', Tau(end));

%% Simulation options

time_current = 0;
t_max = 10;
time_step_pos = 1/500;

%% Controller Options
kv  = [10, 0; 0, 10];
kp  = [100, 0; 0, 100];

%% Trajectory Generation
t = linspace(time_current, t_max, t_max/time_step_pos);
x0 = [0.2,0]';
xf = [0.8,0.3]';
xd = x0+ (xf-x0)* (6*(t/t_max).^5 -15* (t/t_max).^4 + 10* (t/t_max).^3) ;
vd = (xf-x0)/t_max* (30*(t/t_max).^4 -60* (t/t_max).^3 + 30* (t/t_max).^2) ;
ad = (xf-x0)/t_max^2 * (120*(t/t_max).^3 -180* (t/t_max).^2 + 60* (t/t_max)) ;

figure(2)
plot(t,xd(1,:),'r')
hold on
plot(t,vd(1,:),'k')
hold on
plot(t,ad(1,:),'b')
hold on
legend('xd','vd','ad');title('Jerk minimum trajectory');
grid on ; xlabel('Time (s)') ;

%% Plotting Initialization


%% Plant Simulation
x_0 = [xd(:,1);vd(:,1)];

%% zip the data 
 data.fixed_point = fixed_point;
 data.G = [0; -m*g];
 data.tau_max = tau_max;
 data.tau_min = tau_min;
 data.m = m;
 data.kv = kv;
 data.kp = kp;
 data.x0 = x0;
 data.xf = xf; 
 data.t_max = t_max;
%% plot desired trajectory
t = linspace(time_current, t_max, t_max/time_step_pos);
x0 = [0.2,0]';
xf = [0.8,0]';
xds = x0+ (xf-x0)* (6*(t/t_max).^5 -15* (t/t_max).^4 + 10* (t/t_max).^3) ;
vds = (xf-x0)/t_max* (30*(t/t_max).^4 -60* (t/t_max).^3 + 30* (t/t_max).^2) ;
% ads = (xf-x0)/t_max^2 * (120*(t/t_max).^3 -180* (t/t_max).^2 + 60* (t/t_max)) ;

figure(3)
plot(t,xds(1,:),'r')
hold on
plot(t,vds(1,:),'k')
hold on
%% Solving Dynamical Equations
odeopts = odeset('RelTol',1e-5,'AbsTol',1e-5);
[t, x] = ode45(@odefun, [0 t_max], x_0, odeopts, data);
plot(t,x(:,1))
hold on
plot(t,x(:,3))
hold on  

legend('xd','vd','x','v');title('Dynamic Control of three cables system');
grid on ; xlabel('Time (s)') ;

%% Plot the error
figure(4) 
[xds, vds, ads] = traj(t,data);
plot(t,x(:,1:2) - xds');
hold on
plot(t,x(:,3:4) - vds');
legend('e_x','e_y','e_{vx}','e_{vy}');title('Error of the model');
grid on ; xlabel('Time (s)') ;

%% 3D animation
figure(8)

time_step_dynamic = t_max;

x_dynamic = x(:,1);
y_dynamic = x(:,2);


%plot desired and actual static lines
plot(xds(1,:),xds(2,:),'--k')

hold on
plot(x(:,1),x(:,2),'b')
plot(fixed_point(1,:),fixed_point(2,:), 'ro')
grid on
xlim([0, 1]);
ylim([0, 2]);
set(gca,'Color','none');
set(gca,'CLim',[0, 1E-4]);
xlabel("x")
ylabel("y")

%initialize live plots
x_dr = x_dynamic(1);
y_dr = y_dynamic(1);

plot_num = plot(x_dynamic(1),y_dynamic(1),'ko');
plot_cable1 = plot([x_dynamic(1), fixed_point(1,1)], [y_dynamic(1),fixed_point(2,1)],'-y','LineWidth',0.5);
plot_cable2 = plot([x_dynamic(1), fixed_point(1,2)], [y_dynamic(1),fixed_point(2,2)],'-m','LineWidth',0.5);
plot_cable3 = plot([x_dynamic(1), fixed_point(1,3)], [y_dynamic(1),fixed_point(2,3)],'-g','LineWidth',0.5);

for idx_dynamic = 1:length(x_dynamic) 
    

    %plot load
    plot_num.XData = x_dynamic(idx_dynamic); 
    plot_num.YData = y_dynamic(idx_dynamic); 
    
    % plot cable
    plot_cable1.XData = [x_dynamic(idx_dynamic), fixed_point(1,1)]; 
    plot_cable1.YData = [y_dynamic(idx_dynamic), fixed_point(2,1)]; 
    plot_cable2.XData = [x_dynamic(idx_dynamic), fixed_point(1,2)]; 
    plot_cable2.YData = [y_dynamic(idx_dynamic), fixed_point(2,2)]; 
    plot_cable3.XData = [x_dynamic(idx_dynamic), fixed_point(1,3)]; 
    plot_cable3.YData = [y_dynamic(idx_dynamic), fixed_point(2,3)]; 
    pause(0.002)
     
end

%% odefun
function dx = odefun(t,x,data)

    %% extract data
    fixed_point = data.fixed_point;
    G = data.G;
    tau_max = data.tau_max;
    tau_min = data.tau_min;
    xL = x(1:2);
    vL = x(3:4);
    m = data.m;
    kv = data.kv;
    kp = data.kp;
    %% tension distribution based on desired value Tau_ff
    [xLd,vLd,aLd] = traj(t,data);
    Wd = normc(fixed_point - xLd*ones(1, max(size(fixed_point))));
    Nd = null(Wd);
    x0d = Wd\-G;
    lb = tau_min*ones(max(size(fixed_point)),1) - x0d;
    ub = tau_max*ones(max(size(fixed_point)),1) - x0d;
    fun = @(x) norm(x0d+ Nd*x);
    cond_A = [Nd; -Nd];
    cond_b = [ub; -lb];
    lamda = fmincon(fun,-1,cond_A,cond_b);
    Tau_ff = x0d+ Nd*lamda;
    
    %% calculate the controller part Tau_c
    W = normc(fixed_point - xL*ones(1, max(size(fixed_point))));
    Tau_c = W\(m* (aLd + kv* (vLd - vL) + kp* (xLd - xL))); 
    Tau_mg = Tau_ff+ Tau_c;
    
    % Saturation maintain the minimum tension along the cable
    
    if(Tau_mg(1)<1)
        Tau_mg(1) = 1;
    end
    if(Tau_mg(2)<1)
        Tau_mg(2) = 1;
    end
    if(Tau_mg(3)<1)
        Tau_mg(3) = 1;
    end
    
    if(Tau_mg(1)>30)
        Tau_mg(1) = 1;
    end
    if(Tau_mg(2)>30)
        Tau_mg(2) = 1;
    end
    if(Tau_mg(3)>30)
        Tau_mg(3) = 1;
    end
    Tau_mg
    
    
    %% Dynamics 
    vL_dot = (W* Tau_mg + G)/m ;
    xL_dot = vL ;
    dx = [xL_dot; vL_dot] ;
end

%% define a minimum jerk trajectory from [0.2;0] to [0.8;0]
function [xd, vd, ad] = traj(t,data) 
    t_max = data.t_max;
    x0 = data.x0';
    xf = data.xf';
    xd = x0 + (6*(t/t_max).^5 -15* (t/t_max).^4 + 10* (t/t_max).^3)*(xf-x0);
    xd = xd';
    vd = 1/t_max* (30*(t/t_max).^4 -60* (t/t_max).^3 + 30* (t/t_max).^2)*(xf-x0);
    vd = vd';
    ad = 1/t_max^2 * (120*(t/t_max).^3 -180* (t/t_max).^2 + 60* (t/t_max))*(xf-x0);
    ad = ad';
end