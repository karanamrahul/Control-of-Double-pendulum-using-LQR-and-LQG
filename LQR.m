clc
clear

disp("Question D")

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

disp("The eigen value of the controllability matrix is:")

eigs(double(A))

disp("The rank of the controllability matrix is:")
rank(Ctrl)

lamda = eye(6,6) * sym('lambda');


A_lbd = [lamda - A B];

rank(A_lbd);

intial_state = [5;0;10;0;20;0]; 

Q=[10 0 0 0 0 0;
   0 5 0 0 0 0;
   0 0 1000 0 0 0;
   0 0 0 1000 0 0;
   0 0 0 0 100 0;
   0 0 0 0 0 10];
R = 0.0001;

C = eye(6);
D = 0;

sys_rep = ss(A,B,C,D);

K_val = lqr(A,B,Q,R);

sys_rep_k = ss(A-B*K_val,B,C,D);

eigs(A-B*K_val)

figure
initial(sys_rep_k,intial_state)