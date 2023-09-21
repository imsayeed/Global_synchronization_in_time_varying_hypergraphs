
#include <stdio.h>  
#include <stdlib.h> 
#include <math.h> 
#include <time.h> 
int main(void)
{
srand(time(NULL));
int i,j,j1,ii,n,m,mm,j2,j3,i2,i3,i4,i5,l1,rf,rp,np,k1,IRAND,TT,k2,k3,k4,k5,nn,tn,i1,i7,na,d1,d2,N,n1;
float cum[100],v[100],gsc[100],p,temp,t0,tp0,h,l,w,sigma,rho,beta,eps,eta,RAND,ev1,alpha,a,b,c;
float ev2,znorm[100],mle,maxlyap,x10,v20,y10,n20,z10,s20,x1,v2,y1,n2,z1,s2;
float ax10,ay10,az10,ax1,ay1,az1,eps1,eps2,gc2,sig,ax20,ax2,ay20,ay2,az20,az2,ax30,ax3,ay30,ay3,az30,az3,omega;
float k14,k15,k16,k24,k25,k26,k34,k35,k36,k44,k45,k46,k42,k43,k11,k12,k13,k21,k22,k23,k31,k32,k33,k41,gc;
float tau,taus,gca,gk,gs,eca,ek,es,thetas,einh,sigr,sigc,betar,betac;
float k17,k18,k19,k27,k28,k29,k37,k38,k39,k47,k48,k49;
float sig1c,sig1r,sig2c,sig2r,ev3,alpha2,t;

FILE *fopen(),*fp1; ///,*fp2;

fp1=fopen("omega=2_ev2_ev3_sf1.out","w");
///fp2=fopen("L_ts.out","w");
//N=100;
n=2; 
n1=2*n; //
nn=n1*(n1+1); //n*(n+1);
TT=1; 
eps1=1.0;     
alpha2=2.0;
//ev2=1.0;ev3=2.0; 
sigr=1.0;sigc=4.3;betar=1.0;betac=1.1;

omega=2.0;
w=sigc-betac*sigr/betar;
//ev2=1.0;ev3=2.0;
eps1=1.0;
mm=500000; // Time Iteration
tn=200000; // Transient Time

for(i2=0;i2<=200;i2++)
{
ev2=pow(10.0,(-2.0+3.0/200*i2));

for(i3=0;i3<=200;i3++)
{
ev3=pow(10.0,(-2.0+3.0/200*i3));


 mle=0.0;
printf("%d\t%d\n",i2,i3);

for(i1=1;i1<=TT;i1++)
{
//x10= (float)rand()/RAND_MAX; //-1.50+3.50*(float)rand()/RAND_MAX ;/// -0.07+0.05*(float)rand()/RAND_MAX; ///(float)rand()/RAND_MAX ;/// (-1.5+rand()*3.50)*0.01
//y10= (float)rand()/RAND_MAX; //-7.0+8.0*(float)rand()/RAND_MAX ;///0.12*(float)rand()/RAND_MAX; (-7.0+rand()*8.0)*0.01
//z10= (float)rand()/RAND_MAX; //2.9+0.5*(float)rand()/RAND_MAX ;///0.165+0.035*(float)rand()/RAND_MAX;  (2.9+rand()*0.5)*0.01

for(i=n1+1;i<=nn;i++)
v[i]=0.0;

for(i=1;i<=n1;i++)
{
v[(n1+1)*i]=1.0;
cum[i]=0.0;
}

h=0.01; t0=0.0;
for(ii=1;ii<=mm;ii++)
{

//l= sig*(y10-x10);
//k11=l*h;
//l=x10*(rho-z10)-y10;
//k12=l*h;
 
//x1=x10+k11/2.0; 
//y1=y10+k12/2.0; 

//l= sig*(y1-x1);
//k21=l*h;
//l=x1*(rho-z1)-y1;
//k22=l*h;
 
//x1=x10+k21/2.0; 
//y1=y10+k22/2.0; 


//l= sig*(y1-x1);
//k31=l*h;
//l=x1*(rho-z1)-y1;
//k32=l*h;

//x1=x10+k31; 
//y1=y10+k32; 


//l= sig*(y1-x1);
//k41=l*h;
//l=x1*(rho-z1)-y1;
//k42=l*h;

  
//x10=x10+(k11+2.0*(k21+k31)+k41)/6.0;
//y10=y10+(k12+2.0*(k22+k32)+k42)/6.0;
t0=t0+h;

if(ii>tn)
{  //// fprintf(fp2,"%f\t%f\t%f\t%f\t\n",w0,x10,y10,z10);

for(j1=0;j1<=n1-1;j1++)
{

ax20=v[n1+1+j1]; 
ay20=v[2*n1+1+j1]; 
ax30=v[3*n1+1+j1]; 
ay30=v[4*n1+1+j1]; 



l= omega*ax30-2*sigr*ax20-eps1*3*(sigr/betar)*ev2*(cos(2*w*t0)*ax20-sin(2*w*t0)*ay20); 
k14=l*h;
l=omega*ay30-2*betac*(sigr/betar)*ax20-eps1*3*(sigr/betar)*ev2*(sin(2*w*t0)*ax20+cos(2*w*t0)*ay20); 
 k15=l*h;


l= -omega*ax20-2*sigr*ax30-eps1*3*(sigr/betar)*ev3*(cos(2*w*t0)*ax30-sin(2*w*t0)*ay30); 
k17=l*h;
l=-omega*ay20-2*betac*(sigr/betar)*ax30-eps1*3*(sigr/betar)*ev3*(sin(2*w*t0)*ax30+cos(2*w*t0)*ay30); 
 k18=l*h;


ax2=ax20+k14/2.0; 
ay2=ay20+k15/2.0; 

ax3=ax30+k17/2.0; 
ay3=ay30+k18/2.0; 

t=t0 +h/2.0;


l= omega*ax3-2*sigr*ax2-eps1*3*(sigr/betar)*ev2*(cos(2*w*t)*ax2-sin(2*w*t)*ay2); 
k24=l*h;
l=omega*ay3-2*betac*(sigr/betar)*ax2-eps1*3*(sigr/betar)*ev2*(sin(2*w*t)*ax2+cos(2*w*t)*ay2); 
 k25=l*h;


l= -omega*ax2-2*sigr*ax3-eps1*3*(sigr/betar)*ev3*(cos(2*w*t)*ax3-sin(2*w*t)*ay3); 
k27=l*h;
l=-omega*ay2-2*betac*(sigr/betar)*ax3-eps1*3*(sigr/betar)*ev3*(sin(2*w*t)*ax3+cos(2*w*t)*ay3); 
 k28=l*h;

ax2=ax20+k24/2.0;
ay2=ay20+k25/2.0;

ax3=ax30+k27/2.0;
ay3=ay30+k28/2.0;

t=t0+h/2.0;


l= omega*ax3-2*sigr*ax2-eps1*3*(sigr/betar)*ev2*(cos(2*w*t)*ax2-sin(2*w*t)*ay2); 
k34=l*h;
l=omega*ay3-2*betac*(sigr/betar)*ax2-eps1*3*(sigr/betar)*ev2*(sin(2*w*t)*ax2+cos(2*w*t)*ay2); 
 k35=l*h;


l= -omega*ax2-2*sigr*ax3-eps1*3*(sigr/betar)*ev3*(cos(2*w*t)*ax3-sin(2*w*t)*ay3); 
k37=l*h;
l=-omega*ay2-2*betac*(sigr/betar)*ax3-eps1*3*(sigr/betar)*ev3*(sin(2*w*t)*ax3+cos(2*w*t)*ay3); 
 k38=l*h;

ax2=ax20+k34;
ay2=ay20+k35;

ax3=ax30+k37;
ay3=ay30+k38;

t=t0+h;


l= omega*ax3-2*sigr*ax2-eps1*3*(sigr/betar)*ev2*(cos(2*w*t)*ax2-sin(2*w*t)*ay2); 
k44=l*h;
l=omega*ay3-2*betac*(sigr/betar)*ax2-eps1*3*(sigr/betar)*ev2*(sin(2*w*t)*ax2+cos(2*w*t)*ay2); 
 k45=l*h;


l= -omega*ax2-2*sigr*ax3-eps1*3*(sigr/betar)*ev3*(cos(2*w*t)*ax3-sin(2*w*t)*ay3); 
k47=l*h;
l=-omega*ay2-2*betac*(sigr/betar)*ax3-eps1*3*(sigr/betar)*ev3*(sin(2*w*t)*ax3+cos(2*w*t)*ay3); 
 k48=l*h;


ax20=ax20+(k14+2.0*(k24+k34)+k44)/6.0;
ay20=ay20+(k15+2.0*(k25+k35)+k45)/6.0;

ax30=ax30+(k17+2.0*(k27+k37)+k47)/6.0;
ay30=ay30+(k18+2.0*(k28+k38)+k48)/6.0;

t0=t0+h;

v[n1+1+j1]=ax20; 
v[2*n1+1+j1]=ay20; 
v[3*n1+1+j1]=ax30; 
v[4*n1+1+j1]=ay30; 

}
znorm[1]=0.0;
for(j=1;j<=n1;j++)
znorm[1]=znorm[1]+pow(v[n1*j+1],2.0);
znorm[1]=sqrt(znorm[1]);
for(j=1;j<=n1;j++)
v[n1*j+1]=v[n1*j+1]/znorm[1];
for(j=2;j<=n1;j++)
{
for(k1=1;k1<=j-1;k1++)
{
gsc[k1]=0.0;
for(l1=1;l1<=n1;l1++)
gsc[k1]=gsc[k1]+v[n1*l1+j]*v[n1*l1+k1];
}
for(k2=1;k2<=n1;k2++)
{
for(l1=1;l1<=j-1;l1++)
v[n1*k2+j]=v[n1*k2+j]-gsc[l1]*v[n1*k2+l1];
}
znorm[j]=0.0;
for(k3=1;k3<=n1;k3++)
znorm[j]=znorm[j]+pow(v[n1*k3+j],2.0);
znorm[j]=sqrt(znorm[j]);
for(k4=1;k4<=n1;k4++)
v[n1*k4+j]=v[n1*k4+j]/znorm[j];
}
for(k5=1;k5<=n1;k5++)
cum[k5]=cum[k5]+log(znorm[k5])/log(2.0);
}
t0=t0+h;
}// Time Iteration (ii)
tp0=t0-tn*h;
for(i=1;i<=nn;i++)
{
for(j=i+1;j<=nn;j++)
{
if(cum[j]>cum[i])
{
temp=cum[i]; cum[i]=cum[j]; cum[j]=temp;
}
}
}
mle=mle+cum[1]/tp0;//printf("%f\n",cum[1]/t0);
}
fprintf(fp1,"%f\t%f\t%f\n",ev2,ev3,mle/TT);
///printf("%f\t%f\t%f\t%f\n",eps,eps2,eta,mle/TT);
}
}

return 0;
fclose(fp1); //// fclose(fp2);

}


