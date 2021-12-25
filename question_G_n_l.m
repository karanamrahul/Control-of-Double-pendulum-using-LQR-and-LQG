clear
clc

syms M m1 m2 l1 l2 g;

A=[0 1 0 0 0 0; 
    0 0 -(m1*g)/M 0 -(m2*g)/M 0;
    0 0 0 1 0 0;
    0 0 -((M+m1)*g)/(M*l1) 0 -(m2*g)/(M*l1) 0;
    0 0 0 0 0 1;
    0 0 -(m1*g)/(M*l2) 0 -(g*(M+m2))/(M*l2) 0];

B=[0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];

A = double(subs(A, {M, m1, m2, l1, l2, g}, {1000, 100, 100, 20, 10, 9.8}));

B = double(subs(B, {M, m1, m2, l1, l2, g}, {1000, 100, 100, 20, 10, 9.8}));

M=1000;
m1=100;
m2=100;
l1=20;
l2=10;
g=9.81;

Q=[10 0 0 0 0 0;
   0 100 0 0 0 0;
   0 0 10 0 0 0;
   0 0 0 100 0 0;
   0 0 0 0 10 0;
   0 0 0 0 0 10000];
R = 0.01;

initial_x = [0; 0; pi/6; 0; pi/4; 0; 0; 0; 0; 0; 0; 0];
tspan = 0:0.1:5000;

[final_t, final_x] = ode45(@(t, x)cart_system(x, M, m1, m2, l1, l2, g, A, B, Q, R), tspan, initial_x);
hold on
plot(final_t, final_x)
% ylim([-5, 5])
xlabel("time");
ylabel("States")
legend('x', 'v', 'theta1', 'angular1', 'theta2', 'angular2')
hold off

function dx = cart_system(x, M, m1, m2, l1, l2, g, A, B, Q, R)

C1 = [1 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0]; 

[K, ~, ~] = lqr(A,B,Q,R);
F=-K*x(1:6);

Vd=0.1*eye(6);
Vn=5;
Kp = lqr(A',C1',Vd,Vn)';

S =(A-Kp*C1)*x(7:12);

dx=zeros(12,1);
dx(1) = x(2); 
dx(2)=(F-(g/2)*(m1*sind(2*x(3))+m2*sind(2*x(5)))-(m1*l1*(x(4)^2)*sind(x(3)))-(m2*l2*(x(6)^2)*sind(x(5))))/(M+m1*((sind(x(3)))^2)+m2*((sind(x(5)))^2));
dx(3)= x(4);
dx(4)= (dx(2)*cosd(x(3))-g*(sind(x(3))))/l1'; 
dx(5)= x(6); 
dx(6)= (dx(2)*cosd(x(5))-g*(sind(x(5))))/l2; 
dx(7)= x(2)-x(10); 
dx(8)= dx(2)-S(2);
dx(9)= x(4)-x(11);
dx(10)= dx(4)-S(4);
dx(11)= x(6)-x(12);
dx(12)= dx(6)-S(6);
end