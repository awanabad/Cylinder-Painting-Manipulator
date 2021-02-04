function T_mat = T_matrix(q)

L=[1.5, 1.75, 0.3]; %Enter Link Lengths 

%DH Parameters
a = [0, L(2), L(3)];
d = [L(1), 0, 0];
alpha = [-pi/2, 0, 0]; 
theta = [q(1), q(2), q(3)];

%DH Matrix 
A= {1:4, rand(4,4)}; %memory preallocation

for i= 1:numel(a)
    A{i}= [ cos(theta(i)), -sin(theta(i)) * cos(alpha(i)), sin(theta(i))*sin(alpha(i)), a(i)*cos(theta(i)); 
            sin(theta(i)), cos(theta(i))*cos(alpha(i)), -cos(theta(i))*sin(alpha(i)), a(i)*sin(theta(i));
            0, sin(alpha(i)), cos(alpha(i)), d(i); 
            0, 0, 0, 1; ]; 
end

%T Matrices
T_mat={ A{1};  A{1}*A{2};  A{1}*A{2}*A{3} };

end
