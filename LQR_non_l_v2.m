clear
clc

M=1000;
m1=100;
m2=100;
l1=20;
l2=10;
g=9.8;

A=[0 1 0 0 0 0; 
    0 0 -(m1*g)/M 0 -(m2*g)/M 0;
    0 0 0 1 0 0;
    0 0 -((M+m1)*g)/(M*l1) 0 -(m2*g)/(M*l1) 0;
    0 0 0 0 0 1;
    0 0 -(m1*g)/(M*l2) 0 -(g*(M+m2))/(M*l2) 0];

B=[0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];

% Q=[10 0 0 0 0 0;
%    0 100 0 0 0 0;
%    0 0 0 10 0 0;
%    0 0 0 100 0 0;
%    0 0 0 0 10 0;
%    0 0 0 0 0 10000];
% R = 0.01;
Q=[1 0 0 0 0 0;
   0 1 0 0 0 0;
   0 0 1 0 0 0;
   0 0 0 1 0 0;
   0 0 0 0 1 0;
   0 0 0 0 0 1];
R = 1;

initial_x = [2; 0; pi/6; 0; pi/4; 0];
tspan = 0:0.1:5000;
rank(ctrb(A,B))
[K_val, ~, ~] = lqr(A,B,Q,R);

F=@(x)-K_val*x;
eigs(A-B*K_val)

[final_t, final_x] = ode45(@(t, x)cart_system(x, M, m1, m2, l1, l2, g, F(x)), tspan, initial_x);
plot(final_t, final_x)
ylim([-5, 5])
xlabel("time");
ylabel("States")
legend('x', 'v', 'theta1', 'angular1', 'theta2', 'angular2')

function dx = cart_system(x, M, m1, m2, l1, l2, g, F)



dx=zeros(6,1);
dx(1) = x(2); 
dx(2)=(F-(g/2)*(m1*sind(2*x(3))+m2*sind(2*x(5)))-(m1*l1*(x(4)^2)*sind(x(3)))-(m2*l2*(x(6)^2)*sind(x(5))))/(M+m1*((sind(x(3)))^2)+m2*((sind(x(5)))^2));
dx(3)= x(4);
dx(4)= (dx(2)*cosd(x(3))-g*(sind(x(3))))/l1'; 
dx(5)= x(6); 
dx(6)= (dx(2)*cosd(x(5))-g*(sind(x(5))))/l2; 
end