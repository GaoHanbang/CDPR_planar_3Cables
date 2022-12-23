clear
clc
shift_x = 0.2;
shift_y = 0.2;
ts = 2;
syms x(t) y(t) pos(t)
x(t) = piecewise(in(t, 'real') & t ~= 3, t^5/2160 - t^4/144 + t^3/36);
y(t) = piecewise(t < 3, (4*t^5)/3645 - (5*t^4)/486 + (2*t^3)/81, 3 < t, (2*t^3)/81 - (256*(t/6 - 1/2)^5)/15 - (5*t^4)/486 + (4*t^5)/3645);

dx(t) = diff(x,t);
dy(t) = diff(y,t);


t = linspace(0,6,600);
plot(t,x(t));
hold on
plot(t,y(t));

fun_x = @(t) (x(t-ts)+shift_x);
fun_y = @(t) (y(t-ts)+shift_y);
% t = linspace(2,8,600);
% plot(fun_x(t), fun_y(t))

% pos(t) = [fun_x(t);fun_y(t)];
% pos(4)
% t = linspace(2,8,600);
% plot(t, pos(t))
% axis equal

data.x0 = [0.2;0];
data.xf = [0.2;0.2];
data.t_max = 2;
t = linspace(0,2,600);
% [xd, vd, ad] = traj(t',data);
% [xd, vd, ad] = traj(1,data)
% plot(t, traj(t',data))
function x = traj(t,data) 
    t_max = data.t_max;
    x0 = data.x0';
    xf = data.xf';
    xd = x0 + (6*(t/t_max).^5 -15* (t/t_max).^4 + 10* (t/t_max).^3)*(xf-x0);
    xd = xd';
    vd = 1/t_max* (30*(t/t_max).^4 -60* (t/t_max).^3 + 30* (t/t_max).^2)*(xf-x0);
    vd = vd';
    x= [xd;vd];
end