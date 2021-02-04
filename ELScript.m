%************************************************************************
%Euler-Langrange Method
%************************************************************************
clc 
clear all
%declaring manipulator parameters 
syms q1 q2 q3 L1 L2 L3 Lc1 Lc2 Lc3 m1 m2 m3 g %symbols for various parameters 
syms I1xx I1yy I1zz I2xx I2yy I2zz I3xx I3yy I3zz %inertia matrix symbols 
syms qd1 qd2 qd3 qdd1 qdd2 qdd3  %joint velocities and accelerations
q= [q1, q2, q3].'; %joint angles
qd= [qd1, qd2, qd3].'; %joint velocities
qdd= [qdd1, qdd2, qdd3].'; %joint accelerations 
m= [m1, m2, m3]; %mass per link (g) 
L= [L1, L2, L3]; %link lengths (m) 
L_c=[Lc1, Lc2, Lc3]; %distances from link start to centre of mass (m) 
%intertia matrices about centre of masses
I1=[ I1xx, 0, 0; 0, I1yy, 0; 0, 0, I1zz];
I2=[ I2xx, 0, 0; 0, I2yy, 0; 0, 0, I2zz];
I3=[ I3xx, 0, 0; 0, I3yy, 0; 0, 0, I3zz];
I={I1, I2, I3}; 

%************************************************************************
%Computing DH Matrices
a = [0, L(2), L(3)]; %a values 
d = [L(1), 0, 0]; %d values 
a_c = [0, L_c(2), L_c(3)]; %a values, modified for centre of mass
d_c = [L_c(1), 0, 0]; %d values, modified for centre of mass
alpha = [-pi/2, 0, 0]; %alpha angles
theta = [q(1), q(2), q(3)]; %joint angles

%DH Matrix formulation 
A= {1:3, rand(4,4)}; %memory preallocation
A_c= {1:3, rand(4,4)}; %memory preallocation
for i= 1:numel(q) %for loop that computes both A matrices based on link lengths and mass centres
    
    A{i}= [ cos(theta(i)), -sin(theta(i)) * cosd(alpha(i)*(180/pi)), sin(theta(i))*sind(alpha(i)*(180/pi)), a(i)*cos(theta(i)); 
            sin(theta(i)), cos(theta(i))*cosd(alpha(i)*(180/pi)), -cos(theta(i))*sind(alpha(i)*(180/pi)), a(i)*sin(theta(i));
            0, sind(alpha(i)*(180/pi)), cosd(alpha(i)*(180/pi)), d(i); 
            0, 0, 0, 1; ]; 
        
    A_c{i}= [ cos(theta(i)), -sin(theta(i)) * cosd(alpha(i)*(180/pi)), sin(theta(i))*sind(alpha(i)*(180/pi)), a_c(i)*cos(theta(i)); 
            sin(theta(i)), cos(theta(i))*cosd(alpha(i)*(180/pi)), -cos(theta(i))*sind(alpha(i)*(180/pi)), a_c(i)*sin(theta(i));
            0, sind(alpha(i)*(180/pi)), cosd(alpha(i)*(180/pi)), d_c(i); 
            0, 0, 0, 1; ]; 
    
end

%DH matrices 
A1=A{1};
A2=A{1}*A{2};
A3=A{1}*A{2}*A{3};

%DH matrices, modified for centre of mass
Ac1=A_c{1};
Ac2=A{1}*A_c{2};
Ac3=A{1}*A{2}*A_c{3};
Ac={Ac1, Ac2, Ac3};
%************************************************************************
%Computing modified Jacobians
%z vectors
z0=[0 0 1].';
z1= A1(1:3, 3);
z2= A2(1:3, 3);
z3= A3(1:3, 3);
%origin vectors
o0=[0 0 0].';
o1= A1(1:3, 4);
o2= A2(1:3, 4);
o3= A3(1:3, 4);
%centre of mass origin vectors 
oc1=Ac1(1:3, 4);
oc2=Ac2(1:3, 4);
oc3=Ac3(1:3, 4);
%Jacobians
J1=[ cross(z0,oc1 - o0), zeros(3,1), zeros(3,1);
     z0, zeros(3,1), zeros(3,1) ];
J2=[ cross(z0,oc2 - o0), cross(z1,oc2 - o1), zeros(3,1);
    z0, z1, zeros(3,1) ];
J3=[ cross(z0,oc3 - o0), cross(z1,oc3 - o1), cross(z2,oc3 - o2);
    z0, z1, z2 ];
J = {J1, J2, J3};
%************************************************************************
%Computing KEs from Jacobians, Computing Intertia Matrix
Kl={1:3,zeros(3,3)}; %memory preallocation
KW={1:3,zeros(3,3)}; %memory preallocation
d_temp={1:3,zeros(3,3)}; %memory preallocation

for j=1:numel(q)
Jv_temp=J{j}(1:3,1:3); 
Jw_temp=J{j}(4:6,1:3); 
R_temp=Ac{j}(1:3,1:3);
I_temp=I(j);

Kl{j}=m(j)*(Jv_temp).'*Jv_temp;
KW{j}=(Jw_temp).'*R_temp*I_temp*(R_temp).'*Jw_temp;

d_temp{j}=(Kl{j})+(KW{j});
end

D=simplify(d_temp{1}+d_temp{2}+d_temp{3});
%************************************************************************
%Computing Christoffels 
for page=1:numel(q)
    for row=1:numel(q)  
        for col=1:numel(q) 
         
            i_temp=col; %temporary notations to fix swap
            j_temp=row;
            k_temp=page;
            
            d_kj=D(k_temp,j_temp); %retrieving d values 
            d_ki=D(k_temp,i_temp);
            d_ij=D(i_temp,j_temp);
            
            term_1=0.5*diff(d_kj,q(i_temp)); %computing partial derivatives 
            term_2=0.5*diff(d_ki,q(j_temp));
            term_3=0.5*diff(d_ij,q(k_temp));
         
            C(i_temp, j_temp, k_temp)=term_1+term_2+term_3; %summing for C value
            
        end 
    end
    
end

%************************************************************************
%Computing PE
for k=1:numel(q)
    
    PE(k)=m(k)*g*Ac{k}(3,4); 
    G(k)=diff(PE(k),q(k));

end
%************************************************************************
%Finally, Computing Dynamic Equations of Motion
Tao_D=D*qdd;
Tao_C=simplify((C(:,:,1)*qd)+(C(:,:,2)*qd)+(C(:,:,3)*qd).'*qd);
Tao_G=G.';
Tao=simplify(Tao_D+Tao_C+Tao_G);

%************************************************************************






