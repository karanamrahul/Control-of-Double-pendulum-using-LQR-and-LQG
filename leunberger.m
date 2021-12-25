clc
clear

disp("Question F")

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

% Ctrl = ctrb(A, B)
Ctrl = ([B A*B (A^2)*B (A^3)*B (A^4)*B (A^5)*B]);

C1 = [1 0 0 0 0 0; 0 0 0 0 0 0;0 0 0 0 0 0];  
C2 = [0 0 1 0 0 0; 0 0 0 0 1 0;0 0 0 0 0 0]; 
C3 = [1 0 0 0 0 0; 0 0 0 0 1 0;0 0 0 0 0 0]; 
C4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];

D = 0;

Q=[10000 0 0 0 0 0;
   0 1000 0 0 0 0;
   0 0 100 0 0 0;
   0 0 0 10 0 0;
   0 0 0 0 1000 0;
   0 0 0 0 0 100];
R = 0.001;

new_poles=[-0.1;-0.3;-0.5;-0.7;-0.9;-1.1];

K=lqr(A,B,Q,R);

x_initial = [0, 0, pi/4, 0, pi/3, 0, 0, 0, 0, 0, 0, 0];

Luenberger_B = [B;zeros(size(B))];

L1 = place(A',C1',new_poles)'; 
Luenberger_A1 = [(A-B*K) B*K; 
        zeros(size(A)) (A-L1*C1)];
Luenberger_C1 = [C1 zeros(size(C1))];
sys_rep1 = ss(Luenberger_A1,Luenberger_B,Luenberger_C1,D);
figure 
initial(sys_rep1,x_initial)
figure
step(sys_rep1)

L3 = place(A',C3',new_poles)'; 
Luenberger_A3 = [(A-B*K) B*K; 
        zeros(size(A)) (A-L3*C3)];
Luenberger_C3 = [C3 zeros(size(C3))];
sys_rep3 = ss(Luenberger_A3,Luenberger_B,Luenberger_C3,D);
figure 
initial(sys_rep3,x_initial)
figure
step(sys_rep3)

L4 = place(A',C4',new_poles)'; 
Luenberger_A4 = [(A-B*K) B*K; 
        zeros(size(A)) (A-L4*C4)];
Luenberger_C4 = [C4 zeros(size(C4))];
sys_rep4 = ss(Luenberger_A4,Luenberger_B,Luenberger_C4,D);
figure 
initial(sys_rep4,x_initial)
figure
step(sys_rep4)