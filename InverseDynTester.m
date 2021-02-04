%InvDyn Tester

%Test Inputs 
q_test=[0,0,0].';
qd_test=[0,0,0].';
qdd_test=[0,0,0].';

%Function Call
Tao_test=InverseDynamicsCalculator(q_test,qd_test,qdd_test)