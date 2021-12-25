clc
clear

disp("Question E")

syms M m1 m2 l1 l2 g;

A=[0 1 0 0 0 0; 
    0 0 -(m1*g)/M 0 -(m2*g)/M 0;
    0 0 0 1 0 0;
    0 0 -((M+m1)*g)/(M*l1) 0 -(m2*g)/(M*l1) 0;
    0 0 0 0 0 1;
    0 0 -(m1*g)/(M*l2) 0 -(g*(M+m2))/(M*l2) 0];

B=[0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];


C1 = [1 0 0 0 0 0];  
C2 = [0 0 1 0 0 0; 0 0 0 0 1 0]; 
C3 = [1 0 0 0 0 0; 0 0 0 0 1 0]; 
C4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];

O1 = [C1; C1*A; C1*A^2; C1*A^3; C1*A^4; C1*A^5];
disp("The rank for the first C")
rank(O1)

O2 = [C2; C2*A; C2*A^2; C2*A^3; C2*A^4; C2*A^5];
disp("The rank for the second C")
rank(O2)

O3 = [C3; C3*A; C3*A^2; C3*A^3; C3*A^4; C3*A^5];
disp("The rank for the third C")
rank(O3)

O4 = [C4; C4*A; C4*A^2; C4*A^3; C4*A^4; C4*A^5];
disp("The rank for the fourth C")
rank(O4)