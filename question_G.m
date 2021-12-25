clc
clear

disp("Question G")

syms M m1 m2 l1 l2 g;

A=[0 1 0 0 0 0; 
   0 0 -(m1*g)/M 0 -(m2*g)/M 0;
   0 0 0 1 0 0;
   0 0 -((M+m1)*g)/(M*l1) 0 -(m2*g)/(M*l1) 0;
   0 0 0 0 0 1;
   0 0 -(m1*g)/(M*l2) 0 -(g*(M+m2))/(M*l2) 0];

A = double(subs(A, {M, m1, m2, l1, l2, g}, {1000, 100, 100, 20, 10, 9.8}));

B=[0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];

B = double(subs(B, {M, m1, m2, l1, l2, g}, {1000, 100, 100, 20, 10, 9.8}));

C1 = [1 0 0 0 0 0];  
C2 = [0 0 1 0 0 0; 0 0 0 0 1 0]; 
C3 = [1 0 0 0 0 0; 0 0 0 0 1 0]; 
C4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];

D = 0;

Q=[10000 0 0 0 0 0;
   0 1000 0 0 0 0;
   0 0 100 0 0 0;
   0 0 0 10 0 0;
   0 0 0 0 1000 0;
   0 0 0 0 0 100];
R = 0.05;

new_poles=[-0.1;-0.3;-0.5;-0.7;-0.9;-1.1];

K=lqr(A,B,Q,R);

x_initial = [0, 0, pi/4, 0, pi/3, 0, 0, 0, 0, 0, 0, 0];

Vd = 0.1*eye*(6);
Vn = 0.5;

K_C1 = lqr(A', C1', Vd, Vn)';
sys_rep1 = ss([(A-B*K) B*K; zeros(size(A)) (A-K_C1*C1)], [B;zeros(size(B))],[C1 zeros(size(C1))], D);
figure 
initial(sys_rep1,x_initial)
figure
step(sys_rep1)

K_C3 = lqr(A', C3', Vd, Vn)';
sys_rep3 = ss([(A-B*K) B*K; zeros(size(A)) (A-K_C3*C3)], [B;zeros(size(B))],[C3 zeros(size(C3))], D);
figure 
initial(sys_rep3,x_initial)
figure
step(sys_rep3)

K_C4 = lqr(A', C4', Vd, Vn)';
sys_rep4 = ss([(A-B*K) B*K; zeros(size(A)) (A-K_C4*C4)], [B;zeros(size(B))],[C4 zeros(size(C4))], D);
figure 
initial(sys_rep4,x_initial)
figure
step(sys_rep4)