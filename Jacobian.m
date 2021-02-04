function F = Jacobian(q)
%Enter Link Lengths 
L_1=1.5;
L_2=1.75;
L_3=0.3;

%Enter Theta Values
%DH Parameters
a = [0, L_2, L_3];
d = [L_1, 0, 0];
alpha = [-pi/2, 0, 0]; 
theta = [q(1), q(2), q(3)];

%Creating Homogenous Transformation Matrices 
for i= 1:numel(a)

    A{i}= [ cos(theta(i)), -sin(theta(i)) * cos(alpha(i)), sin(theta(i))*sin(alpha(i)), a(i)*cos(theta(i)); 
            sin(theta(i)), cos(theta(i))*cos(alpha(i)), -cos(theta(i))*sin(alpha(i)), a(i)*sin(theta(i));
            0, sin(alpha(i)), cos(alpha(i)), d(i); 
            0, 0, 0, 1; ] ;
end
%Creating Tranformation Matrix
T02= A{1} * A{2};
T03= A{1} * A{2} * A{3}; 
%Creating Jacobian Matrix 
 o3 = T03(1:3, 4);
 z0 = [0; 0; 1];
 o0 = [0; 0; 0];
 J1 = [cross(z0, (o3-o0)); z0];
 A1=A{1};
 z1 = A1(1:3, 3);
 o1 = A1(1:3, 4);
 J2 = [cross(z1, (o3-o1)); z1];
 z2 = T02(1:3, 3);
 o2 = T02(1:3, 4);
 J3 = [cross(z2, (o3-o2)); z2];
 J = [J1, J2, J3];
 
 %Output Data
disp('***************************************************');
disp('The Jacobian Matrix is');
disp('***************************************************');  
disp('J='); 
disp(J);
disp('***************************************************'); 
end