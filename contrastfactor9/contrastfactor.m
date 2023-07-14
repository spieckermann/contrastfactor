#! /usr/bin/octave -qf
#Calculation of Contrast Factor for arbitrary symetry materials
#See Martinez-Garcia 2009,Acta Cryst A65 pp 109-119 and related
#Version 9.0.3
#2023/07/14

#Published under CC by 4.0
#
#Usage: Edit Material in the appropriate file 'material.m' then run as:-> contrastfactor.m material
#
#scalars are in  lower case, vectors and tensors upper case

# please refer to https://doi.org/10.25365/thesis.12708
#and to  DOI: 10.5281/zenodo.8146313 ( latest DOI  https://zenodo.org/badge/latestdoi/666299173)

#set as script
1;

%use gnuplot
graphics_toolkit gnuplot 

%Functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plane
% calculate "x" or "y" values for the normal vector "normvec"
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = planey(x,z,normvec)
	y=-(normvec(1)*x+normvec(3)*z)/normvec(2);
endfunction
function x = planex(y,z,normvec)
	x=-(normvec(2)*y+normvec(3)*z)/normvec(1);
endfunction
function x = planez(x,y,normvec)
	x=-(normvec(1)*x+normvec(2)*y)/normvec(3);
endfunction



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Preamble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
printf('\n----------------------------------\n');
printf('- Calculation of Contrast Factor -\n');
printf('- Florian Spieckermann 2009-2011 -\n');
printf('----------------------------------\n\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Deal with command line arguments 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arg_list= argv ();

if (nargin==1)
	#load Material Constants
	material=arg_list{1};
	[Stiffness,a,b,c,alpha,beta,gamma,hkl_matrix,uvw_matrix,disltype_matrix,diffr_vector_matrix]=eval(material);
else
	printf('The contrastfactor calculation requires one argumnet');
endif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate Relevant Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%orthonormal system
I_orth=[1;0;0];
J_orth=[0;1;0];
K_orth=[0;0;1];

% Volume of Unit cell
V=a*b*c*sqrt(1-cosd(alpha)^2-cosd(beta)^2-cosd(gamma)^2+2*cosd(alpha)*cosd(beta)*cosd(gamma)); %V=det(inv(M));

%reciprocal parameters in ortho coords
a_star=b*c*sind(alpha)/V;
b_star=c*a*sind(beta)/V;
c_star=a*b*sind(gamma)/V;
alpha_star=asin(V/(a*b*c*sind(beta)*sind(gamma)))/pi*180;
beta_star=asin(V/(a*b*c*sind(alpha)*sind(gamma)))/pi*180;
gamma_star=asin(V/(a*b*c*sind(alpha)*sind(beta)))/pi*180;
fflush(stdout);

%definition of M
M=[1/a,0,0;
   -cosd(gamma)/(a*sind(gamma)),1/(b*sind(gamma)),0;
   a_star*cosd(beta_star),b_star*cosd(alpha_star),c_star];

%calculation of G_m from M
G_m=inv(M)*inv(M');

%Vector formulation of unit cell in orthonormal variables
A=inv(M)'*I_orth;
B=inv(M)'*J_orth;
C=inv(M)'*K_orth;

A_star=cross(B,C)/V;
B_star=cross(C,A)/V;
C_star=cross(A,B)/V;

%Calculate Twotheta for plot
lambda=1.542; %wavelength (Angstrom)
hd=diffr_vector_matrix(:,1);
kd=diffr_vector_matrix(:,2);
ld=diffr_vector_matrix(:,3);
printf('Diffraction Vector:\n');
printf('reflection \t twotheta\n');
%general formulation for twoteta values
twotheta=zeros(rows(hd),1);
for dnr=1:rows(hd)
twotheta(dnr) =2*asind(lambda/2*sqrt((hd(dnr).*A_star+kd(dnr).*B_star+ld(dnr).*C_star)'*(hd(dnr).*A_star+kd(dnr).*B_star+ld(dnr).*C_star)));
printf('(%g %g %g) \t %8.4f \n',hd(dnr),kd(dnr),ld(dnr),twotheta(dnr));
endfor
%twotheta
fflush(stdout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cycle over slip systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialise Contrastfactor CF and angle between Burgers vector and diffraction vector Delta
CF=[];
Delta=[];
for hklnr=1:rows(hkl_matrix)
printf('*******************************************************************\n');
printf('* Slip System: slip plane: (%g %g %g) burgers vector:[%g %g %g], phi= %g.\n',hkl_matrix(hklnr,:),uvw_matrix(hklnr,:),disltype_matrix(hklnr));
printf('*******************************************************************\n');

%slip plane (reciprocal coordinates)
HKL=[hkl_matrix(hklnr,1);hkl_matrix(hklnr,2);hkl_matrix(hklnr,3)];
h=HKL(1);
k=HKL(2);
l=HKL(3);
%burgersvector (lattice coordinates)
UVW=[uvw_matrix(hklnr,1);uvw_matrix(hklnr,2);uvw_matrix(hklnr,3)];
u=UVW(1);
v=UVW(2);
w=UVW(3);
fflush(stdout);
%dislocation character (0 screw, 90 edge)
phi=disltype_matrix(hklnr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculation of transformation matrix P=[C_1,C_2,C_3]
%(unit vectors of slip system C_1, C_2, C_3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%C_2
%slip plane normal N (orto coords)
N=h*A_star+k*B_star+l*C_star;
%absolute length of N
n=norm(N);
%coordinates of C_2
C_2=1/n*M*HKL;

%C_3
% Burgers Vector Burg (ortho coords)
Burg=u*A+v*B+w*C;
%absolute length burgers vector
burg=sqrt(UVW'*G_m*UVW);
% normalized burgers vector in ortho coords
Cb_3=Burg/burg;
%define rotation matrix Rot(phi,C_2)
c21=C_2(1);
c22=C_2(2);
c23=C_2(3);
Rot=[c21^2*(1-cosd(phi))+cosd(phi),c21*c22*(1-cosd(phi))+c23*sind(phi),c21*c23*(1-cosd(phi))-c22*sind(phi);
c21*c22*(1-cosd(phi))-c23*sind(phi),c22^2*(1-cosd(phi))+cosd(phi),c22*c23*(1-cosd(phi))+c21*sind(phi);
c21*c23*(1-cosd(phi))+c22*sind(phi),c22*c23*(1-cosd(phi))-c21*sind(phi),c23^2*(1-cosd(phi))+cosd(phi)];

C_3=Rot*Cb_3;

%C_1
C_1=cross(C_2,C_3);

%Definition of transformation matrix P
P=[C_1';C_2';C_3'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate E_ijmn through eigenvalue problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%transform Stiffness to slip system coords
TStiffness=voigt2tensor(Stiffness); 	%voigt 2 tensor notation
TCturn=transform(TStiffness,(inv(P))');	%transformation ortho to slip coordinate system Stiffness  (P transforms from ortho to slip)
Cturn=tensor2voigt(TCturn);		%transform back to voigt
% define Q, R ,T see Ting aniostropic elasticity
Q=reshape(TCturn(:,1,:,1),3,3);
R=reshape(TCturn(:,1,:,2),3,3);
T=reshape(TCturn(:,2,:,2),3,3);
%Define Nfund fundamental Elasticity matrix (see Ting)
N1=-inv(T)*R';
N2=inv(T);
N3=R*inv(T)*R'-Q;
N4=N1';
Nfund=[N1,N2;N3,N4];
%save Nfund of slip system
%fname=["Nfund",num2str(hklnr),".mat"]
%save("-mat",fname,"Nfund")
% check if eigenvalue problem is ill posed
%conditionnumber=cond(Nfund)
fflush(stdout);

% solve eigenvalue problem with eig
[EVECT,EVAL]=eig(Nfund);
%normalize eigenvectors
EVAL = EVAL(:,1) ./ sum(EVAL(:,1));
psol=eig(Nfund);
A_alpha=EVECT(1:3,:);
L_alpha=EVECT(4:6,:);
%norm(A_alpha(:,1))
fflush(stdout);


Burgnorm=(inv(M)*inv(P))'*UVW/burg;

%calculate E_ijmn
phi_polar=vec(0:0.0001:2*pi(1,"double"));	
%initialize Etensor
Etensor=ones(3,2,3,2);
countere=1;
for eindi=1:3
	for eindj=1:2
		for eindm=1:3
			for eindn=1:2
				countere++; %setting all 36 values of matrix?
				%normalized burgersvector in slip system coordinates : (inv((P*M)')*UVW)/burg #check!!!
				beta_ij = dstrain(A_alpha,L_alpha,psol,Burgnorm,phi_polar,eindi,eindj);
				beta_mn = dstrain(A_alpha,L_alpha,psol,Burgnorm,phi_polar,eindm,eindn);
				#plot(phi_polar,cumtrapz(phi_polar,beta_ij.*beta_mn))
				%pause(0.01);
				Evalu=1/pi*trapz(phi_polar,beta_ij.*beta_mn);
				Etensor(eindi,eindj,eindm,eindn)=Evalu;
			endfor
		endfor
	endfor	
endfor

% transform E_ijmn to 6*6 matrix E_hat
E_hat=ones(6,6);
countereh=1;
for tri=1:3
	for trj=1:2
		for trm=1:3
			for trn=1:2
				countereh++;
				if trj==1
					rk=tri;
				endif
				if trj==2
					rk=tri+3;
				endif
				if trn==1
					rl=trm;
				endif
				
				if trn==2
	  			rl=trm+3;
				endif
				E_hat(rk,rl)=Etensor(tri,trj,trm,trn);
				endfor
		endfor
	endfor	
endfor
printf('E_hat=\n');
printf('| %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f |\n',E_hat)
fflush(stdout);

%Nfund*[A_alpha(:,1);L_alpha(:,1)]-psol(1)*[A_alpha(:,1);L_alpha(:,1)]


%calculate burgers normal plane for plotting
%if Burg(2)==0
%	zv1 = [ -c : 1 : c ]; yv1 = [ -b : 1 : b ];
%	[zzv1,yyv1] = meshgrid(zv1,yv1);
%	xxv1 = planex(yyv1,zzv1,Burg);
%	%printf('case1\n')
%elseif (Burg(1)==0) || (Burg(3)==0)
%	zv1 = [ -a : 1 : c ]; xv1 = [ -a : 1 : a ];
%	[zzv1,xxv1] = meshgrid(zv1,xv1);
%	yyv1 = planey(xxv1,zzv1,Burg);
%	%printf('case2\n')
%else
%	xv1 = [ -a : 1 : a ]; yv1 = [ -b : 1 : b ];
%	[xxv1,yyv1] = meshgrid(xv1,yv1);
%	zzv1 = planez(xxv1,yyv1,Burg);
%	%printf('case3\n')
%endif
%mesh(xxv1,yyv1,zzv1)
%hold on



%plotting coordinate systems %orthonormal is given by frame
%plot3([0;A(1)],[0;A(2)],[0;A(3)],"2","linewidth",2,[0;B(1)],[0;B(2)],[0;B(3)],"2","linewidth",3,[0;C(1)],[0;C(2)],[0;C(3)],"2;Crystal;","linewidth",2,[0;C_1(1)],[0;C_1(2)],[0;C_1(3)],"3","linewidth",2,[0;C_2(1)],[0;C_2(2)],[0;C_2(3)],"3;Slip;","linewidth",2,[0;C_3(1)],[0;C_3(2)],[0;C_3(3)],"1;C3;","linewidth",2,[0;Burg(1)],[0;Burg(2)],[0;Burg(3)],"4;Burg;","linewidth",2);
%pltitle=sprintf('# Slip System: (%g %g %g) [%g %g %g] phi= %g\n Reflex: (%g %g %g) -> chi=%4.1f °, eta=%4.1f °',HKL,UVW,phi,  Diffr,chi, eta );
%title(pltitle);
%xlabel("I_{ortho}")
%ylabel("J_{ortho}")
%zlabel("K_{ortho}")
%axis([-12,12,-12,12,-12,12])
%hold off
%pause(0.1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Start CF Calculation for the different reflexes in the slipsystem 
%(cycle over diffr_vector matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CF_sys=[];
Delta_sys=[];
for diffnr=1:rows(diffr_vector_matrix)

% calculate diffraction vector in different coordinate systems
Diffr=[diffr_vector_matrix(diffnr,1);diffr_vector_matrix(diffnr,2);diffr_vector_matrix(diffnr,3)];  %in the reciprocal frame
DIFFR=Diffr(1)*A_star+Diffr(2)*B_star+Diffr(3)*C_star; %in orthonormal coordinates 
%define absloute length of diffracrion vector
diffr=norm(DIFFR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate Geo_ijmn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau1=DIFFR'*C_1/(diffr*norm(C_1));
tau2=DIFFR'*C_2/(diffr*norm(C_2));
tau3=DIFFR'*C_3/(diffr*norm(C_3));
tau=[tau1;tau2;tau3];

%calculate Geo_ijmn
%initialize Geo
Geo=ones(3,2,3,2); 
counterg=1;
for gindi=1:3
	for gindj=1:2
		for gindm=1:3
			for gindn=1:2
				counterg++; %setting all 36 values of matrix?
				gval=tau(gindi)*tau(gindj)*tau(gindm)*tau(gindn);
				Geo(gindi,gindj,gindm,gindn)=gval;
				fflush (stdout);
			endfor
		endfor
	endfor	
endfor

% transform Geo_ijmn to 6*6 matrix Geo_hat
Geo_hat=ones(6,6);
countergh=1;
for tri=1:3
	for trj=1:2
		for trm=1:3
			for trn=1:2
				countergh++;
				if trj==1
					rk=tri;
				endif
				if trj==2
					rk=tri+3;
				endif
				if trn==1
					rl=trm;
				endif
				
				if trn==2
	  			rl=trm+3;
				endif
				Geo_hat(rk,rl)=Geo(tri,trj,trm,trn);
				endfor
		endfor
	endfor	
endfor
%printf('Geo_hat=\n');
%printf('| %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f |\n',Geo_hat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate Contrastfactor(hkl) from G_ijmn E_ijmn and (see e.g. martinez-garcia2009)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initailize
cfactor=0;
counterc=1;
for conti=1:3
	for contj=1:2
		for contm=1:3
			for contn=1:2
			counterc++;
			etval=Etensor(conti,contj,contm,contn);
			gtval=Geo(conti,contj,contm,contn);
			contrval=Etensor(conti,contj,contm,contn)*Geo(conti,contj,contm,contn);
			cfactor=cfactor+contrval;
			fflush(stdout);
			endfor
		endfor
	endfor	
endfor
CF_sys=[CF_sys;cfactor];

%Calculate relevant angles delta,chi, eta

%%direction cosines between Diffr and C-Axis
delta=acosd(dot(DIFFR,C)/(norm(DIFFR)*norm(C)));
%%direction cosines between Diffr and B-Axis
%delta=acosd(dot(DIFFR,B)/(norm(DIFFR)*norm(B)));
%%direction cosines between Diffr and A-Axis
%delta=acosd(dot(DIFFR,A)/(norm(DIFFR)*norm(A)));

Delta_sys=[Delta_sys;delta];


%%direction cosines between Diffr and Line(chi) and Burgers Vector (eta)
chi=atan2(norm(cross(DIFFR,C_3)),dot(DIFFR,C_3))*180/pi();
%eta=atan2(norm(cross(DIFFR,Burg)),dot(DIFFR,Burg))*180/pi();
eta=acosd(dot(DIFFR,Burg)/(norm(DIFFR)*burg));

% Print Result
printf('The Contrastfactor of the (%g %g %g) Reflex is: %8.4f (chi=%4.1f°, eta=%4.1f°) \n',  Diffr , cfactor, chi,eta);
%end diffr_vector cycling
endfor

% store CF as matrix for later usage
CF=[CF,CF_sys];
Delta=[Delta,Delta_sys];

%end slip system cycling
endfor

printf('\n');
fflush(stdout);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate results for all slip systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

system("mkdir results");

%Save CF 
%octave format 
savefile_oct=['./results/CF_octave_',material,'.dat'];
save ("-text",savefile_oct,"diffr_vector_matrix","twotheta","hkl_matrix","uvw_matrix","disltype_matrix","CF")
% nice table for gnumeric or excel
[status, username]=system("echo -n $(whoami)");
[status, hostname]=system("echo -n $(hostname)");
extime=strftime ('%r (%Z) %A %e %B %Y', localtime (time ()));
savefile=['./results/CF_table_',material,'.dat'];
printf('Saving Results to %s.\n',savefile);
fd = fopen(savefile,'w');
fprintf(fd,'# %s for %s produced by %s@%s %s\n',savefile, material,username,hostname,extime);
fprintf(fd,'#Peakindex\t2theta(°)\t');
for ind=1:rows(hkl_matrix) 
fprintf(fd,'C_[%g%g%g](%g%g%g)_%g\t',hkl_matrix(ind,:),uvw_matrix(ind,:),disltype_matrix(ind));
endfor
fprintf(fd,'\n');
for ind=1:rows(CF)
	fprintf(fd,'[%2g %2g %2g]\t%7.3g\t', diffr_vector_matrix(ind,:), twotheta(ind));
	for ind2=1:columns(CF) 
	fprintf(fd,'%10.8g\t',CF(ind,ind2));
	endfor
	fprintf(fd,'\n');
endfor
fclose(fd);

%plot CF(delta)
for ind=1:columns(CF)
	subplot(columns(CF),1,ind)
	[dsort_sys,idsort_sys]=sort(Delta(:,ind));
	plot(dsort_sys,CF(:,ind)(idsort_sys),"-o")
	titlestring=sprintf('Slip System (%g %g %g) [%g %g %g], phi= %g',hkl_matrix(ind,:),uvw_matrix(ind,:),disltype_matrix(ind));
	title(titlestring);
	xlabel("{/Symbol d} (degree)");
	ylabel("CF")
endfor
plotfile_delta=['./results/Delta_',material,'.eps'];
print(plotfile_delta,"-depsc2","-F24")
pause
subplot(1,1,1)

