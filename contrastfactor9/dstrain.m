%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dstrain
% calculates element nm of derivative of strain field beta_mn: name: dstrain at polar variable phi_polar
% check if al matrix operations are really correct
% Burgnorm ......... normalized burgers vector in slip system coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beta_mn = dstrain(A_alpha,L_alpha,psol,Burgnorm,phi_polar,m,n)
phi_polar=vec(phi_polar);
beta_alp=[];
for alpha=1:2:6 %using alpha values 1,3,5 because imaginary values must be positive (TING page 137)
	D_alpha=-(L_alpha(:,alpha).'*Burgnorm)/(A_alpha(:,alpha).'* L_alpha(:,alpha)); %for the dot product we have to use the unitary (complex) space inner product definition 
	beta_val=(A_alpha(m,alpha)*D_alpha*psol(alpha)^(n-1))./(cos(phi_polar)+psol(alpha)*sin(phi_polar));
	beta_alp=[beta_alp,beta_val];
endfor

%plot(phi_polar,imag(beta_alp(:,1)),phi_polar,imag(beta_alp(:,2)),phi_polar,imag(beta_alp(:,3)));
%pause(0.1)

beta_mn=imag(beta_alp(:,1)+beta_alp(:,2)+beta_alp(:,3));

%set pind=1 to plot phi_polar dependence and pind=2 to plot 3D (0 for no plot)
pind=2;

if (pind==1)
%plot
plot(phi_polar,beta_mn,"-1",phi_polar,zeros(rows(phi_polar),1)*beta_mn(1))
beta_title=sprintf('{/Symbol b}_{%g %g}\n',m,n);
beta_ylabel=sprintf('{/Symbol b}_{%g %g}\n',m,n);
%title(beta_title);
xlabel('angle {/Symbol f} (radians)');
ylabel(beta_ylabel);
axis([0 2*pi()])
betafilename=sprintf('./results/plots/beta%g%g.pdf',m,n);
print(betafilename,"-dpdfcairo","-FHelvetics:28")
pause(0.1)


elseif (pind==2)
%3dplot

x=-10:0.5:10;
y=-10:0.5:10;
[xx,yy]=meshgrid(x,y);
beta_alp2=0;
[phi_polar2,r]=cart2pol(xx,yy);
for alpha=1:2:6
	D_alpha2=-L_alpha(:,alpha).'*Burgnorm/(A_alpha(:,alpha).'*L_alpha(:,alpha));
	beta_val2=A_alpha(m,alpha)*D_alpha2*psol(alpha)^(n-1)./(cos(phi_polar2)+psol(alpha)*sin(phi_polar2));
	beta_alp2=beta_alp2 + beta_val2;
endfor

#surf(x,y,imag(beta_alp2));
surf(x,y,imag(beta_alp2),'FaceColor','interp',...
   'EdgeColor',[0.4 0.4 0.4],...
   'FaceLighting','gouraud')
%axis tight
%view(-50,30)
%camlight left


%contourf(x,y,imag(beta_alp2))
%view(90,90)
%axis x "square"
%axis y "square"
beta_title=sprintf('{/Symbol b}_{%g %g}\n',m,n);
betafilename=sprintf('./results/plots/beta%g%g_3D.pdf',m,n)
beta_zlabel=sprintf('{/Symbol b}_{%g %g}\n',m,n);

xlabel('x (Å)');
ylabel('y (Å)');
zlabel(beta_zlabel);
%title(beta_title);
print(betafilename,"-dpdfcairo","-FHelvetics:22")
pause(0.01)

endif



endfunction

