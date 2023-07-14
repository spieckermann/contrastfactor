%Cu
function [Stiffness,a,b,c,alpha,beta,gamma,hkl_matrix,uvw_matrix,disltype_matrix,diffr_vector_matrix]=Cu

% Material Parameters!
Stiffness=[16.84,12.14,12.14,0,0,0;
12.14,16.84,12.14,0,0,0;
12.14,12.14,16.84,0,0,0;
0,0,0,7.54,0,0;
0,0,0,0,7.54,0;
0,0,0,0,0,7.54]
fflush(stdout);
a=3.6148
b=a
c=a
alpha=90
beta=90
gamma=90
fflush(stdout);

%Dislocation Parameters
%hkl_matrix=[0,1,0;
%1,0,0;
%0,1,0;
hkl_matrix=[0,1,0;
1,0,0;
0,1,0] 
fflush(stdout);
%uvw_matrix=[0,0,1;
%0,0,1; 
%1,0,0;
uvw_matrix=[0,0,1;
0,0,1; 
1,0,0]
fflush(stdout);

%dislocation character (0 screw, 90 edge
%disltype_matrix=[0;
%0;
%0;
disltype_matrix=[90;
90;
90]
fflush(stdout);


%Diffraction Vector
diffr_vector_matrix=[1,1,1;
2,0,0;
2,2,0;
3,1,1;
2,2,2;
4,0,0];
fflush(stdout);





endfunction
