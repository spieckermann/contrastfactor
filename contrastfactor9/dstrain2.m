%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dstrain
% calculates element nm of derivative of strain field beta_mn: name: dstrain at polar variable phi_polar
% check if all matrix operations are really correct
% Burgnorm ......... normalized burgers vector in slip system coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beta_mn = dstrain2(A_alpha,L_alpha,psol,Burgnorm,phi_polar,m,n)
phi_polar=vec(phi_polar);
beta_alp=[];
dot_arg_conj_prefer=2;
for alpha =1:2:6 %using alpha values 1,3,5 because imaginary values must be positive (positive definite solutions TING page 137)
	%test=dot(L_alpha(:,alpha),Burgnorm)-dot(Burgnorm,L_alpha(:,alpha))
	D_alpha=-(L_alpha(:,alpha).'*Burgnorm)/(A_alpha(:,alpha).'* L_alpha(:,alpha)); %for the dot product we have to use the unitary (complex) space inner product definition 
	%value in the sum
	beta_val=(A_alpha(m,alpha)*D_alpha*(psol(alpha)).^(n-1))./(cos(phi_polar)+psol(alpha)*sin(phi_polar));
	beta_alp=[beta_alp,beta_val];
endfor
beta_mn=imag(beta_alp(:,1)+beta_alp(:,2)+beta_alp(:,3));
%plot(phi_polar,beta_mn)
%pause(0.01)

endfunction
