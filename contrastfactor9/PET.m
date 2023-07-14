%PET not complete!!!!
function [Stiffness,a,b,c,alpha,beta,gamma,hkl_matrix,uvw_matrix,disltype_matrix,diffr_vector_matrix]=PET

% Material Parameters PET (Calculated for 300 K Rutledge 1997)
Stiffness=[ 14.4, 6.4, 3.4, -2.2, -0.3, -1.8;
	     6.4,17.3, 9.5,  3.3, -0.5,  0.5; 
             3.4, 9.5, 178,  3.8, -0.7, -1.8;
            -2.2, 3.3, 3.8,  6.6,  0.2, -0.4; 
            -0.3,-0.5,-0.7,  0.2,  1.4,    0;
            -1.8, 0.5,-1.8, -0.4,    0,  1.2]
%lattice parameters PET triclinic (Angstrom)
a=4.56
b=5.94
c=10.75
alpha=98.5
beta=118
gamma=112

%Dislocation Parameters (Ahzi 1994)
hkl_matrix=[1,0,0;
0,1,0;
1,0,0;
1,0,0;
0,1,0;
1,0,0] 
fflush(stdout);

uvw_matrix=[0,0,1;
0,0,1; 
0,1,0;
0,0,1;
0,0,1; 
0,1,0]
fflush(stdout);

%dislocation character (0 screw, 90 edge
disltype_matrix=[0;
0;
0;
90;
90;
90]
fflush(stdout);


%PET
#diffr_vector_matrix=[0,-1,1;
#		     0,1,0;
#                     1,-1,-1;
#                     1,-1,0;
#                     0,1,-3;
#                     1,0,0]
#                     1,0,1]
%PET gerald
diffr_vector_matrix=[0,1,-1;
   0,1,0;
   1,-1,-1;
   1,-1,0;
   1,-1,-2;
   1,0,0;
   1,-1,1]
fflush(stdout);



endfunction
