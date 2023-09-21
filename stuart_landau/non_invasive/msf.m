clear all; clc;

sigr=1.0;sigc=4.3;betar=1.0;betac=1.1;
pp=2*betac*(sigr/betar);
omega=6.0;
x=[];
for i=1:100
    ev2=10^(-2.0+3.0/100*i);
    for j=1:100
        ev3=10^(-2.0+3.0/100*j);
A=[-2*sigr+ev2,0,omega,0;-pp,-ev2,0,omega;-omega,0,-2*sigr+ev3,0;0,-omega,-pp,-ev3];
e=eig(A);
e=real(e);
emax=max(e);
x=[x;ev2,ev3,emax];
    end
end

n=100;%number of oscillators
m=100;% number of itterations
u=zeros(n,m);
%v=u;
v=zeros(n,m);
z=zeros(n,m);
%z1=zeros(n,m);x
for i=1:n
    %u(i,j)=x(i,x1); 
    for j=1:m
        
      u(i,j)=x(m*(i-1)+j,1);
      % if(i==1) v(j)=x(100*(i-1)+j,2);endclc
      v(i,j)=x(m*(i-1)+j,2);  
      %z(i,j)=smooth(x(m*(i-1)+j,3),500);
      z(i,j)=x(m*(i-1)+j,3);
    end
end

%my_map=load('my_map.txt');
%load('MyColormaps','mycmap')
%c=z;
figure()
%cmap=colormap((inferno));  
%subplot(1,2,2)

surface(u,v,z,'Marker','.')
shading interp
axis tight %([0 0.1 0 0.25])
%colorbar;
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
datarange = [-0.01 0.01];
breakpoint = 0.0;
ctlen = 256;
colors = [0 0 0; 0.9714    0.9855    0.6028];
pos = (max(datarange)-breakpoint)/diff(datarange);
n = [floor(ctlen*(1-pos)) ceil(ctlen*pos)];
cmap = [repmat(colors(1,:),[n(1),1]); repmat(colors(2,:),[n(2),1])];
colormap(inferno);
%caxis(datarange)
colorbar
