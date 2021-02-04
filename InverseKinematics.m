function F = InverseKinematics(T)
%Enter Link Lengths 
L_1=1.5;
L_2=1.75;
L_3=0.3;

%Enter Theta Values
syms t1 t2 t3;  
%DH Parameters
a = [0, L_2, L_3];
d = [L_1, 0, 0];
alpha = [-pi/2, 0, 0]; 
theta = [t1, t2, t3];
%Creating Homogenous Transformation Matrices 
for i= 1:numel(a)

    A{i}= [ cos(theta(i)), -sin(theta(i)) * cos(alpha(i)), sin(theta(i))*sin(alpha(i)), a(i)*cos(theta(i)); 
            sin(theta(i)), cos(theta(i))*cos(alpha(i)), -cos(theta(i))*sin(alpha(i)), a(i)*sin(theta(i));
            0, sin(alpha(i)), cos(alpha(i)), d(i); 
            0, 0, 0, 1; ] ;
end
%Tranformation Matrix
Tsym= A{1} * A{2} * A{3}; 

%Equating symbols
eq1 = T(1,4)==Tsym(1,4);
eq2 = T(2,4)==Tsym(2,4);
eq3 = T(3,4)==Tsym(3,4);

[t1,t2,t3] = solve(eq1,eq2,eq3);

disp('***************************************************');
disp('Inverse Kinematics Result');
disp('***************************************************');
disp('Possible theta 1 values are');
disp(t1);
disp('Possible theta 2 values are');
disp(t2);
disp('Possible theta 3 values are');
disp(t3);
disp('***************************************************');
end