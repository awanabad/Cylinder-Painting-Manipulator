%Euler Lagrange Method
%************************************************************************
%Code Section 1: Input Parameters 
%************************************************************************
syms q1 q2 q3 qd1 qd2 qd3 qdd1 qdd2 qdd3 
q= [q1, q2, q3].'; %joint angles
qd= [qd1, qd2, qd3].'; %joint velocities
qdd= [qdd1, qdd2, qdd3].'; %joint accelerations 
m= [945.89, 480.13, 112.70]; %mass per link (g) 
L= [1.5, 1.75, 0.3]; %link lengths (m) 
L_c=[0.65, 0.84, 0.17]; %distances from link start to centre of mass (m)
g=9.81; %grav constant 
%intertia matrices about centre of masses
I1=[ 196.52, 0, 0; 0, 53.6, 0; 0, 0, 191.83 ];
I2=[ 8.41, 0, 0; 0, 144.37, 0; 0, 0, 147.97 ];
I3=[ 3.52, 0, 0; 0, 5.10, 0; 0, 0, 3.15 ];
%************************************************************************
%Code Section 2: Computing DH Matrices 
%************************************************************************
%DH Parameters 
a = [0, L(2), L(3)]; %a values 
d = [L(1), 0, 0]; %d values 
a_c = [0, L_c(2), L_c(3)]; %a values, modified for centre of mass
d_c = [L_c(1), 0, 0]; %d values, modified for centre of mass
alpha = [-pi/2, 0, 0]; %alpha angles
theta = [q(1), q(2), q(3)]; %joint angles
%Regulr DH Matrices  
A1= [ cos(theta(1)), -sin(theta(1)) * cosd(alpha(1)*(180/pi)), sin(theta(1))*sind(alpha(1)*(180/pi)), a(1)*cos(theta(1)); 
            sin(theta(1)), cos(theta(1))*cosd(alpha(1)*(180/pi)), -cos(theta(1))*sind(alpha(1)*(180/pi)), a(1)*sin(theta(1));
            0, sind(alpha(1)*(180/pi)), cosd(alpha(1)*(180/pi)), d(1); 
            0, 0, 0, 1; ]; 
        
Temp_A2= [ cos(theta(2)), -sin(theta(2)) * cosd(alpha(2)*(180/pi)), sin(theta(2))*sind(alpha(2)*(180/pi)), a(2)*cos(theta(2)); 
            sin(theta(2)), cos(theta(2))*cosd(alpha(2)*(180/pi)), -cos(theta(2))*sind(alpha(2)*(180/pi)), a(2)*sin(theta(2));
            0, sind(alpha(2)*(180/pi)), cosd(alpha(2)*(180/pi)), d(2); 
            0, 0, 0, 1; ]; 
Temp_A3= [ cos(theta(3)), -sin(theta(3)) * cosd(alpha(3)*(180/pi)), sin(theta(3))*sind(alpha(3)*(180/pi)), a(3)*cos(theta(3)); 
            sin(theta(3)), cos(theta(3))*cosd(alpha(3)*(180/pi)), -cos(theta(3))*sind(alpha(3)*(180/pi)), a(3)*sin(theta(3));
            0, sind(alpha(3)*(180/pi)), cosd(alpha(3)*(180/pi)), d(3); 
            0, 0, 0, 1; ];  
        
%Modded DH Matrices  
Ac1= [ cos(theta(1)), -sin(theta(1)) * cosd(alpha(1)*(180/pi)), sin(theta(1))*sind(alpha(1)*(180/pi)), a_c(1)*cos(theta(1)); 
            sin(theta(1)), cos(theta(1))*cosd(alpha(1)*(180/pi)), -cos(theta(1))*sind(alpha(1)*(180/pi)), a_c(1)*sin(theta(1));
            0, sind(alpha(1)*(180/pi)), cosd(alpha(1)*(180/pi)), d_c(1); 
            0, 0, 0, 1; ] ;
             
Temp_Ac2= [ cos(theta(2)), -sin(theta(2)) * cosd(alpha(2)*(180/pi)), sin(theta(2))*sind(alpha(2)*(180/pi)), a_c(2)*cos(theta(2)); 
            sin(theta(2)), cos(theta(2))*cosd(alpha(2)*(180/pi)), -cos(theta(2))*sind(alpha(2)*(180/pi)), a_c(2)*sin(theta(2));
            0, sind(alpha(2)*(180/pi)), cosd(alpha(2)*(180/pi)), d_c(2); 
            0, 0, 0, 1; ] ;        
        
Temp_Ac3= [ cos(theta(3)), -sin(theta(3)) * cosd(alpha(3)*(180/pi)), sin(theta(3))*sind(alpha(3)*(180/pi)), a_c(3)*cos(theta(3)); 
            sin(theta(3)), cos(theta(3))*cosd(alpha(3)*(180/pi)), -cos(theta(3))*sind(alpha(3)*(180/pi)), a_c(3)*sin(theta(3));
            0, sind(alpha(3)*(180/pi)), cosd(alpha(3)*(180/pi)), d_c(3); 
            0, 0, 0, 1; ] ;
%Computing Transformation Matrices
A2=A1*Temp_A2;
A3=A1*Temp_A2*Temp_A3;
Ac2=A1*Temp_Ac2;
Ac3=A1*Temp_A2*Temp_Ac3;
%************************************************************************
%Code Section 3: Computing Jacobian Matrices 
%************************************************************************        
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
%************************************************************************
%Code Section 4: Computing Linear and Rotational Kinetic Energy
%************************************************************************
%Velocity Jacobians 
Jv1=J1(1:3,1:3);
Jv2=J2(1:3,1:3);
Jv3=J3(1:3,1:3);
%Linear KE Calculation
KL1=m(1)*(Jv1.')*Jv1;
KL2=m(2)*(Jv2.')*Jv2;
KL3=m(3)*(Jv3.')*Jv3;
%************************************************************************
%Rotational Jacobians 
Jw1=J1(4:6,1:3);
Jw2=J2(4:6,1:3);
Jw3=J3(4:6,1:3);
%R values 
R1=Ac1(1:3,1:3);
R2=Ac2(1:3,1:3);
R3=Ac3(1:3,1:3);
%Rotational KE Calculation
KW1=(Jw1.')*R1*I1*(R1.')*Jw1;
KW2=(Jw2.')*R2*I2*(R2.')*Jw2;
KW3=(Jw3.')*R3*I3*(R3.')*Jw3;
%************************************************************************
%Code Section 5: Computing Inertia Matrix
%************************************************************************
%Summing Kinetic Energies
SumKE1=KL1*KW1;
SumKE2=KL2*KW2;
SumKE3=KL3*KW3;
%Formulating D Matrix 
D=simplify(SumKE1+SumKE2+SumKE3);
%************************************************************************
%Code Section 6: Computing Christoffel Symbols 
%************************************************************************
%d values 
d11=D(1,1);
d22=D(2,2);
d33=D(3,3);

d12=D(1,2);
d13=D(1,3);       
d21=d12;
d31=d13;

d23=D(2,3);
d32=d23;
%************************************************************************
%differential terms 
d11q1=diff(d11,q(1)); %with respect to q1
d22q1=diff(d22,q(1));
d33q1=diff(d33,q(1));

d12q1=diff(d12,q(1));
d13q1=diff(d13,q(1));
d21q1=d12q1;
d31q1=d13q1;

d23q1=diff(d23,q(1));
d32q1=d23q1;
%************************************************************************
d11q2=diff(d11,q(2)); %with respect to q2
d22q2=diff(d22,q(2));
d33q2=diff(d33,q(2));

d12q2=diff(d12,q(2));
d13q2=diff(d13,q(2));
d21q2=d12q2;
d31q2=d13q2;

d23q2=diff(d23,q(2));
d32q2=d23q2;
%************************************************************************
d11q3=diff(d11,q(3)); %with respect to q3
d22q3=diff(d22,q(3));
d33q3=diff(d33,q(3));

d12q3=diff(d12,q(3));
d13q3=diff(d13,q(3));
d21q3=d12q3;
d31q3=d13q3;

d23q3=diff(d23,q(3));
d32q3=d23q3;
%************************************************************************
%cij1 values 
c111=0.5*(d11q1+d11q1-d11q1);
c221=0.5*(d12q2+d12q2-d22q1);
c331=0.5*(d13q3+d13q3-d33q1);

c121=0.5*(d12q1+d11q2-d11q1);
c131=0.5*(d13q1+d11q3-d11q1);

c211=c121;
c311=c131;

c231=0.5*(d13q2+d12q3-d23q1);
c321=c231;

%cij2 values 
c112=0.5*(d21q1+d21q1-d11q2);
c222=0.5*(d22q2+d22q2-d22q2);
c332=0.5*(d23q3+d23q3-d33q2);

c122=0.5*(d22q1+d21q2-d11q2);
c132=0.5*(d23q1+d21q3-d11q2);

c212=c122;
c312=c132;

c232=0.5*(d23q2+d22q3-d23q2);
c322=c232;

%cij3 values 
c113=0.5*(d31q1+d31q1-d11q3);
c223=0.5*(d32q2+d32q2-d22q3);
c333=0.5*(d33q3+d33q3-d33q3);

c123=0.5*(d32q1+d31q2-d11q3);
c133=0.5*(d33q1+d31q3-d11q3);

c213=c123;
c313=c133;

c233=0.5*(d33q2+d32q3-d23q3);
c323=c233;
%************************************************************************
%C matrices
C1_temp=[c111, c211, c311; c121, c221, c321; c131, c231, c331];
C2_temp=[c112, c212, c312; c122, c222, c322; c132, c232, c332];
C3_temp=[c113, c213, c313; c123, c223, c323; c133, c233, c333];

C1=simplify(C1_temp);
C2=simplify(C2_temp);
C3=simplify(C3_temp);
%************************************************************************
%Code Section 7: Computing Potential Energy Gradient  
%************************************************************************
%Potential Energy Calculation
PE1=m(1)*g*Ac1(3,4);
PE2=m(2)*g*Ac2(3,4);
PE3=m(3)*g*Ac3(3,4);
%g values
g1=diff(PE1,q(1));
g2=diff(PE2,q(2));
g3=diff(PE3,q(3));
%Formulating G matrix
G=[g1 g2 g3].';
%************************************************************************
%Code Section 8: Bringing it all together and Computing Tao  
%************************************************************************        
%Computing Individual Tao Terms 
Term1=simplify((D*qdd));
Term2=simplify(((((C1*qd)+(C2*qd)+(C3*qd)).')*qd));
Term3=G;
%Summation and Simplification of Tao Terms 
TaoTemp=simplify(Term1+Term2+Term3); 

        