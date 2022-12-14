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
t_max = 10;
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
grid on ; xlabel('x(m)') ; ylabel('y(m)') ;

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
time_step_pos = 1/100;

%% Controller Options
kv  = [50, 0; 0, 50];
kp  = [200, 0; 0, 200];

%% Trajectory Generation
interval_2 = cal_via_jerk_traj;
%% Plotting Initialization
s = linspace(0,10,1000);
figure(2);
generated_traj = mytraj(s,interval_2);

plot(s,generated_traj(:,1:2));
yyaxis right
plot(s,generated_traj(:,3:4));

legend('x','y','v_x','v_y');title('Trajectory with respect to time(s)');
grid on ; xlabel('Time (s)') ;
yyaxis left; ylabel('position(m)')
yyaxis right; ylabel('velocity(m/s)')
figure(3);
generated_traj = mytraj(s,interval_2);
plot(generated_traj(:,1),generated_traj(:,2))
title('Trajectory in 2D phase');
axis equal
grid on ; xlabel('position (m)') ; ylabel('position (m)') ;
%% Plant Simulation
x_0 = calc_traj(0,interval_2)';
%% zip the data 
 data.fixed_point = fixed_point;
 data.G = [0; -m*g];
 data.tau_max = tau_max;
 data.tau_min = tau_min;
 data.m = m;
 data.kv = kv;
 data.kp = kp;
 data.interval_2 = interval_2;
%  data.x0 = x0;
%  data.xf = xf; 
%  data.t_max = t_max;
%% plot desired trajectory
t = linspace(time_current, t_max, t_max/time_step_pos);
pos_vel = mytraj(t,interval_2);


figure(4)
yyaxis left
plot(t,pos_vel(:,1),'r')
hold on
plot(t,pos_vel(:,2),'b')
hold on
yyaxis right
plot(t,pos_vel(:,3),'k')
hold on
plot(t,pos_vel(:,4),'g')
hold on


%% Tension info
tension_info = []; % tension_info = [Tau_mg;Tau_ff;Tau_c]
save('tension_info.mat','tension_info');
if isfile('tension_info.mat')
    delete('tension_info.mat');
    save('tension_info.mat','tension_info');
else
     save('tension_info.mat','tension_info');
end
%% Solving Dynamical Equations
odeopts = odeset('RelTol',1e-5,'AbsTol',1e-5);
[t, x] = ode45(@odefun, [0 t_max], x_0, odeopts, data);
yyaxis left
plot(t,x(:,1))
hold on
plot(t,x(:,2))
hold on
yyaxis right
plot(t,x(:,3))
hold on  
plot(t,x(:,4))
hold on
legend('x_d','y_d','v_{dx}','v_{dy}','x','y','v_x','v_y');title('Dynamic Control of three cables system');
grid on ; xlabel('Time (s)');
yyaxis left; ylabel('error of position (m)');
yyaxis right; ylabel('error of velocity (m/s)');

%% Plot the error
figure(5) 
res = mytraj(t,interval_2);
plot(t,x(:,1:2) - res(:,1:2));
hold on
plot(t,x(:,3:4) - res(:,3:4));
legend('e_x','e_y','e_{vx}','e_{vy}');title('Error of the model');
grid on ; xlabel('Time (s)') ;

%% Plot the tension with respect to time

% ATTENTION: the problem of this method to plot cable tension is the time
% is not continous. For one given time, there could be more
% corresponding cable tension. Simply because ode45 function will calculate
% backforwardly the previous time period.
figure(6)
load('tension_info.mat','tension_info')
plot(tension_info(1,:),tension_info(2:end,:));
legend('\tau_{1mg}','\tau_{2mg}','\tau_{3mg}','\tau_{1ff}','\tau_{2ff}','\tau_{3ff}','\tau_{1c}','\tau_{2c}','\tau_{3c}');title('Tensions along three cable');
grid on ; xlabel('Time (s)') ; ylabel('cable tension (N)');

figure(7)  
plot(tension_info(1,:),tension_info(2:4,:));
legend('\tau_{1mg}','\tau_{2mg}','\tau_{3mg}');title('Tensions along three cable');
grid on ; xlabel('Time (s)') ; ylabel('cable tension (N)');



%% 3D animation
figure(8)

time_step_dynamic = t_max;

x_dynamic = x(:,1);
y_dynamic = x(:,2);


%plot desired and actual static lines
plot(generated_traj(:,1),generated_traj(:,2),'--k')

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
    interval_2 = data.interval_2;
    %% tension distribution based on desired value Tau_ff
    res = mytraj(t,interval_2);
    xLd = res(:,1:2)';
    vLd = res(:,3:4)';
%     [xLd,vLd,aLd] = traj(t,data);
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
    Tau_c = W\(m* (kv* (vLd - vL) + kp* (xLd - xL))); 
    Tau_mg = Tau_ff+ Tau_c;
%     Tau_mg = Tau_ff;
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
    
    % Control tension should be less than maximum tension
    if(Tau_mg(1)>30)
        Tau_mg(1) = 1;
    end
    if(Tau_mg(2)>30)
        Tau_mg(2) = 1;
    end
    if(Tau_mg(3)>30)
        Tau_mg(3) = 1;
    end
    %% Store tension infomation
     load('tension_info.mat','tension_info');
     tension_info = [tension_info,[t;Tau_mg;Tau_ff;Tau_c]]; 
     save('tension_info.mat','tension_info');
    %% Dynamics 
    vL_dot = (W* Tau_mg + G)/m ;
    xL_dot = vL ;
    dx = [xL_dot; vL_dot] ;
end

%% define a minimum jerk trajectory consists three intervals:

function pos_vel = mytraj(t,interval_2)

    pos_vel=zeros(max(size(t)),4);   

    for i=1:length(t)       
        pos_vel(i,:)= calc_traj(t(i),interval_2)';  
    end

end

function pos_vel = calc_traj(t,interval_2)
    x0 = [0.2;0];
    x01 = [0.2;0.2];
    x02 = [0.8;0.2];
    xf = [0.8;0];
    t_max = 1;

    %% five viapoints trajectory
%     if(t>=0 && t< 1)       
%         pos_vel = cal_jerk_traj(t,x0,x01,t_max) ;
%     elseif(t>= 1 && t<2)
%         pos_vel = [x01;0;0];
%     elseif(t>=2 && t<8)
%         pos_vel = interval_2(t-2)+ [0.2;0.2;0;0];
%     elseif(t> 8 && t<9)
%         pos_vel = [x02;0;0];
%     else
%         pos_vel = cal_jerk_traj(t-9,x02,xf,t_max) ;
%     end
    
    %% only jerk viapoint curve
        if(t >= 0 && t<2)
            pos_vel = [x0;0;0];
        elseif(t>= 2 && t<8)
            pos_vel = interval_2(t-2)+ [0.2;0;0;0];
        else
            pos_vel = [xf;0;0];
        end
%     fprintf('simulating process: t = %f s \n',t)
end


function res = cal_jerk_traj(t,x0,xf,t_max) 
    xd = x0 + (6*(t/t_max).^5 -15* (t/t_max).^4 + 10* (t/t_max).^3)*(xf-x0);
    vd = 1/t_max* (30*(t/t_max).^4 -60* (t/t_max).^3 + 30* (t/t_max).^2)*(xf-x0);
    res = [xd; vd];
    
end

