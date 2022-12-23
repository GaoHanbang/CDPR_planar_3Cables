%% This script simulate the CDPR system with different frequencies.


%% Parameter Initializations
clc
clear
close all
fixed_point = [0,1; 0.5,1.5;1,1]';
max(size(fixed_point));

m = 1; %kg
g = 9.81 ; % gravity 
G = [0; -m*g];
%% Initial condition
p_0 = [0.5,0]';
tau_min = 1;
tau_max = 30;

%% Simulator options 
t_max = 10;
f_simu = 100; % frequency of the simulator 
f_control = 50; % frequency of the controller
dt = 1/f_simu;
time_step_pos = 1/f_control;
time_current = 0;


%% Trajectory Generation
interval_2 = cal_via_jerk_traj;
all_times = linspace(time_current, t_max, t_max* f_control);
pos_vel = mytraj(all_times,interval_2);
xd_generated = pos_vel(:,1:2);
vd_generated = pos_vel(:,3:4);
x_0 = [xd_generated(1:2,:); vd_generated(1:2,:)];


%% Controller Options
kv  = [50, 0; 0, 50];
kp  = [200, 0; 0, 200];


%% zip the data 
 data.fixed_point = fixed_point;
 data.G = [0; -m*g];
 data.tau_max = tau_max;
 data.tau_min = tau_min;
 data.m = m;
 data.kv = kv;
 data.kp = kp;
 data.interval_2 = interval_2;


%% Plotting Initialization

times = zeros(t_max/dt, 1);
pos = zeros(t_max/dt, 2);
vel = zeros(t_max/dt, 2);
tau_mg = zeros(t_max/dt, 3);
tau_c = zeros(t_max/dt, 3);
tau_ff = zeros(t_max/dt, 3);

%% Plant Simulation

index = 1; %
internal_index = 1;
time_prev = time_current;
odeopts = odeset('RelTol',1e-5,'AbsTol',1e-5);
pos_0 = [0;0];
vel_0 = [0;0];
while (time_current < t_max)
    
    if (t_max - time_current < dt*1e-2)
        break;
    end
    
    if (index==1 || ( (time_current-time_prev) >= time_step_pos-0.00000001 ) ) % runs at f_position
        xLd = xd_generated(internal_index, :)';
        vLd = vd_generated(internal_index, :)';
        % to integrate the error for PID controller
%         if (index > 1) 
%             
%                 temp1 = pos(1:2:index-1,:);
% 
%                 temp2 = xd_generated(1:internal_index-1,:);
% 
%                 e_x = temp1 - temp2;
%                 
%                 integral_x = trapz( all_times(1:internal_index-1), e_x, 1 )';
% 
%             else  
%                 integral_x = [0;0;0];
%        end

        % TDA algorithm
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
        
        % PD Controller 
        W = normc(fixed_point - pos_0 *ones(1, max(size(fixed_point))));
        Tau_c = W\(m* (kv* (vLd - vel_0) + kp* (xLd - pos_0))); 
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
        
        time_prev = time_current;
        internal_index = internal_index + 1;
    end
    data.W =  W;
    data.tau_mg = Tau_mg;
    x_0 = [pos_0; vel_0];
    [~, x_var] = ode45(@odefun, [time_current time_current+dt], x_0, odeopts, data); 
    %% Store and update variables
    
    times(index,1) = time_current;
    pos(index,:) = x_var(max(size(x_var)),1:2)';
    vel(index,:) = x_var(max(size(x_var)),3:4)';
    tau_mg(index,:) = Tau_mg;
    tau_ff(index,:) = Tau_ff;
    tau_c(index,:) = Tau_c;
    %update data
    pos_0 = x_var(max(size(x_var)),1:2)';
    vel_0 = x_var(max(size(x_var)),3:4)';
    %update time
    time_current = time_current + dt;
    index = index + 1;
end



%% Solving Dynamical Equations




%% Plot the error



%% 3D animation

%% odefun
function dx = odefun(t,x,data)

    %% extract data
    fixed_point = data.fixed_point;
    G = data.G;
    xL = x(1:2);
    vL = x(3:4);
    m = data.m;
    kv = data.kv;
    kp = data.kp;
    W = data.W;
    Tau_mg = data.tau_mg;
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

