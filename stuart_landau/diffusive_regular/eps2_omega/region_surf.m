 clear all;
 clc;
%x=load('hr_nnlcl_elc_chm_chim_Rn_1_SI_clrbr_nw_dlta_0.04.dat');
%y=inst_freq1;%file name

%x=load('hh_eps_gc.txt');
% for i=1:10201
% if ((x(i,7)==-1.0)) 
%     x(i,7)=15.0;
% end
% end
x=load('eps1=1.2_eps2_omega_3.out');
% for j=1:6561
%     if(x(j,3)==-1.0)
%         x(j,3)=1e5;
%     end
% end
% for j=1:6561
%     if(x(j,4)==-1.0)
%         x(j,4)=1e5;
%     end
% end
n=101;%number of oscillators
m=101;% number of itterations
u=zeros(n,m);
%v=u;
v=zeros(n,m);
z=zeros(n,m);
%z1=zeros(n,m);x
for i=1:n
    %u(i,j)=x(i,x1); 
    for j=1:m
        
      u(i,j)=x(m*(i-1)+j,3);
      % if(i==1) v(j)=x(100*(i-1)+j,2);endclc
      v(i,j)=x(m*(i-1)+j,2);  
      %z(i,j)=smooth(x(m*(i-1)+j,3),500);
      z(i,j)=x(m*(i-1)+j,4);
    end
end

%my_map=load('my_map.txt');
%load('MyColormaps','mycmap')
%c=z;
figure()
cmap=colormap((inferno));  
%subplot(1,2,2)

surface(u,v,z,'Marker','.')
shading interp
axis tight %([0 0.1 0 0.25])
colorbar;
%set(gca, 'YScale', 'log')
