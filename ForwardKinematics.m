function F = ForwardKinematics(q)
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

%Tranformation Matrix
T= A{1} * A{2} * A{3}; 
Text=['For the angles ', num2str(q(1)),', ',  num2str(q(2)), ', and ',  num2str(q(3)), ', '];

%Output Data
disp('***************************************************');
disp('Forward Kinematics Result');
disp('***************************************************'); 
disp(Text); 
disp('the output transformation matrix is'); 
disp('T='); 
disp(T);
disp('***************************************************'); 

end

