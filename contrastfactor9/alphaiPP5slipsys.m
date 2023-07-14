%alpha-iPP monoclinic 
function [Stiffness,a,b,c,alpha,beta,gamma,hkl_matrix,uvw_matrix,disltype_matrix,diffr_vector_matrix]=alphaiPP

% Material Parameters
Stiffness=[7.78,3.91,3.72,0,0.9,0;
3.91,11.55,3.99,0,-0.36,0; 
3.72,3.99,42.44,0, -0.57, 0;
0,0,0,4.02,0,-0.12; 
0.9,-0.36,-0.57,0,3.1,0;
0,0,0,-0.12,0,2.99]
%lattice parameters iPP (Angstrom)
%6.6500  20.9600   6.5000  90.00  99.33  90.00      Cc
a=6.6500
b=20.9600
c=6.5000
alpha=90
beta=99.33
gamma=90

%Dislocation Parameters
hkl_matrix=[0,1,0;
1,0,0;
1,1,0;
0,1,0;
1,0,0;
0,1,0;
1,0,0;
1,1,0;
0,1,0;
1,0,0] 
fflush(stdout);

uvw_matrix=[0,0,1;
0,0,1; 
0,0,1;
1,0,0;
0,1,0;
0,0,1;
0,0,1; 
0,0,1;
1,0,0;
0,1,0]
fflush(stdout);

%dislocation character (0 screw, 90 edge
disltype_matrix=[0;
0;
0;
0;
0;
90;
90;
90;
90;
90]
fflush(stdout);



%a-iPP
diffr_vector_matrix=[1,1,0;
		     0,4,0;
                     1,3,0;
                     1,1,1;
                     1,3,-1;
                     0,4,1;
		     0,6,0;
		     1,3,1;
		     0,2,1];					
fflush(stdout);



endfunction
