clc
clear
x0 = [0;0];
x01 = [0;0.2];
x02 = [0.8;0.2];
xf = [0.8;0];
t_max = 2;
interval_2 = cal_via_jerk_traj;
function mytraj(t)


function res = cal_jerk_traj(t,x0,xf,t_max) 
    xd = x0 + (6*(t/t_max).^5 -15* (t/t_max).^4 + 10* (t/t_max).^3)*(xf-x0);
    vd = 1/t_max* (30*(t/t_max).^4 -60* (t/t_max).^3 + 30* (t/t_max).^2)*(xf-x0);
    res = [xd; vd];
    
end

% t = linspace(0,6,600);
% xd = traj(t);
% plot(t,xd(1,:));
% hold on
% plot(t,xd(2,:));
% function xd = traj(t)
% x = t.^5/2160 - t.^4/144 + t.^3/36;
% y = ((4*t.^5)/3645 - (5*t.^4)/486 + (2*t.^3)/81).*(t>=0&t<3) +...
%     ((2*t.^3)/81 - (256*(t/6 - 1/2).^5)/15 - (5*t.^4)/486 + (4*t.^5)/3645).*(t>=3&t<=6);
% xd = [x;y];
% 
% end
