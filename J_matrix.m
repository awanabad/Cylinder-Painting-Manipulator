%Jacobian Maker
function J_mat = J_matrix(q)
%T Matrices
T=T_matrix(q);
%z vectors
z_0 = [0 0 1]';
z_1= (T{1}(1:3, 3));
z_2= (T{2}(1:3, 3));
%o vectors
o_0= [0 0 0]';
o_1= T{1}(1:3, 4) ;
o_2= T{2}(1:3, 4); 
o_3= T{3}(1:3, 4); 

%Jacobian
J_1=[cross(z_0, (o_1 - o_0)), zeros(3,1), zeros(3,1)]; 
J_2=[cross(z_0, (o_2 - o_0)), cross(z_1, (o_2 - o_0)), zeros(3,1)]; 
J_3=[cross(z_0, (o_3 - o_0)), cross(z_1, (o_3 - o_1)), cross(z_2, (o_3 - o_2))];

J_mat = {J_1', J_2', J_3'};
end