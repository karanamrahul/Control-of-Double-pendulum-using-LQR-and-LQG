clc
clear

disp("Question C")

syms M m1 m2 l1 l2 g;

A=[0 1 0 0 0 0; 
    0 0 -(m1*g)/M 0 -(m2*g)/M 0;
    0 0 0 1 0 0;
    0 0 -((M+m1)*g)/(M*l1) 0 -(m2*g)/(M*l1) 0;
    0 0 0 0 0 1;
    0 0 -(m1*g)/(M*l2) 0 -(g*(M+m2))/(M*l2) 0];

B=[0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];

% Ctrl = ctrb(A, B)
Ctrl = [B A*B (A^2)*B (A^3)*B (A^4)*B (A^5)*B];


disp("The rank of the controllability matrix is:")
rank(Ctrl)

disp("The determinant of the controllability matrix is:")
disp(det(Ctrl));


lamda = eye(6,6) * sym('lambda');


A_lbd = [lamda - A B];

rank(A_lbd);