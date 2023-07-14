%alphaMg2SiO2 orthorhombic 
function [Stiffness,a,b,c,alpha,beta,gamma,hkl_matrix,uvw_matrix,disltype_matrix,diffr_vector_matrix]=alphaMg2SiO2

% Material Parameters
Stiffness=[328,69,69,0,0,0;
69,200,73,0,0,0;
69,73,235,0,0,0;
0,0,0,66.7,0,0;
0,0,0,0,81.3,0;
0,0,0,0,0,80.9];
a=4.775;
b=10.190;
c=5.978;
alpha=90;
beta=90;
gamma=90;
fflush(stdout);

%Dislocation Parameters
hkl_matrix=[0,1,0;
1,0,0;
0,1,0;
0,1,0;
1,0,0;
0,1,0];
fflush(stdout);

uvw_matrix=[0,0,1;
0,0,1; 
1,0,0;
0,0,1;
0,0,1; 
1,0,0];
fflush(stdout);

%dislocation character (0 screw, 90 edge
disltype_matrix=[0;
0;
0;
90;
90;
90];
fflush(stdout);



%a-Mg2SiO2
diffr_vector_matrix=[0,2,0;
1,1,0;
0,2,1;
1,0,1;
1,1,1;
1,2,0;
1,2,1];
fflush(stdout);



endfunction
