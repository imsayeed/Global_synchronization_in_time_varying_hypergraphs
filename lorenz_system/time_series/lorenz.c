
#include <stdio.h>  
#include <stdlib.h> 
#include <math.h> 
#include <time.h> 
int main(void)
{
srand(time(NULL));
int i,j,j1,ii,n,m,mm,j2,j3,i2,i3,i4,i5,l1,rf,rp,np,k1,IRAND,TT,k2,k3,k4,k5,nn,tn,i1,i7,na,d1,d2,N,n1;
float cum[100],v[100],gsc[100],p,temp,t0,tp0,h,l,w,sigma,rho,beta,eps,eta,RAND,ev1,alpha,a,b,c;
float ev2,znorm[100],mle,maxlyap,x10,x20,y10,y20,x30,y30,x1,x2,y1,y2,x3,y3,z10,z20,z30,z1,z2,z3;
float ax10,ay10,az10,ax1,ay1,az1,gc2,sig,ax20,ax2,ay20,ay2,az20,az2,ax30,ax3,ay30,ay3,az30,az3,omega;
float k14,k15,k16,k24,k25,k26,k34,k35,k36,k44,k45,k46,k42,k43,k11,k12,k13,k21,k22,k23,k31,k32,k33,k41,gc;
float tau,taus,gca,gk,gs,eca,ek,es,thetas,einh,sigr,sigc,betar,betac,eps1,eps2;
float k17,k18,k19,k27,k28,k29,k37,k38,k39,k47,k48,k49,x[4],y[4],z[4];
float sig1c,sig1r,sig2c,sig2r,ev3,alpha2,t,pi,a12,a21,a13,a31,a23,a32,a123;

FILE *fopen(),*fp1,*fp2;

fp1=fopen("eps1=0.5_eps2=0.03_omega=3.0_sig1c_sig2c_2.out","w");
fp2=fopen("eps1=0.5_eps2=0.03_Lorenz_omega=3.0_ts.out","w");
//N=100;
n=2; 
n1=2*n; //
nn=n1*(n1+1); //n*(n+1);
TT=1; 
pi=3.141592653589793;     
alpha2=2.0;
ev2=1.0;ev3=2.0; 
sigr=1.0;sigc=4.3;betar=1.0;betac=1.1;

omega=1.0;
w=sigc-betac*sigr/betar;

sig=10.0;beta=8.0/3.0;rho=28.0;

mm=250000; // Time Iteration
tn=245000; // Transient Time

eps1=0.5;//1.0;
eps2=0.03;//0.02;

for(i1=1;i1<=TT;i1++)
{

x10= (float)rand()/RAND_MAX; //-1.50+3.50*(float)rand()/RAND_MAX ;/// -0.07+0.05*(float)rand()/RAND_MAX; ///(float)rand()/RAND_MAX ;/// (-1.5+rand()*3.50)*0.01
y10= (float)rand()/RAND_MAX; //-7.0+8.0*(float)rand()/RAND_MAX ;///0.12*(float)rand()/RAND_MAX; (-7.0+rand()*8.0)*0.01
z10= (float)rand()/RAND_MAX;
x20= (float)rand()/RAND_MAX; //-1.50+3.50*(float)rand()/RAND_MAX ;/// -0.07+0.05*(float)rand()/RAND_MAX; ///(float)rand()/RAND_MAX ;/// (-1.5+rand()*3.50)*0.01
y20= (float)rand()/RAND_MAX; //-7.0+8.0*(float)rand()/RAND_MAX ;///0.12*(float)rand()/RAND_MAX; (-7.0+rand()*8.0)*0.01
z20= (float)rand()/RAND_MAX;
x30= (float)rand()/RAND_MAX; //-1.50+3.50*(float)rand()/RAND_MAX ;/// -0.07+0.05*(float)rand()/RAND_MAX; ///(float)rand()/RAND_MAX ;/// (-1.5+rand()*3.50)*0.01
y30= (float)rand()/RAND_MAX; //-7.0+8.0*(float)rand()/RAND_MAX ;///0.12*(float)rand()/RAND_MAX; (-7.0+rand()*8.0)*0.01
z30= (float)rand()/RAND_MAX;



h=0.01; t0=0.0;
for(ii=1;ii<=mm;ii++)
{

a12=1/2 - cos(pi/3 + 2*omega*t0)/3;
a13= cos(2*omega*t0)/3 + 1/2;
a23=  1/2 - cos(pi/3 - 2*omega*t0)/3;
a21=a12;
a31=a13;
a32=a23;
a123=1 - (2*cos(pi/3 + 2*omega*t0))/3;

l= sig*(y10-x10)+eps2*(a123*(pow(x20,2.0)*x30-pow(x10,3.0))+a123*(pow(x30,2.0)*x20-pow(x10,3.0)));
k11=l*h;
l=x10*(rho-z10)-y10+eps1*(a12*(y20-y10)+a13*(y30-y10));
k12=l*h;
l=x10*y10-beta*z10;
 k13=l*h;
l= sig*(y20-x20)+eps2*(a123*(pow(x10,2.0)*x30-pow(x20,3.0))+a123*(pow(x30,2.0)*x10-pow(x20,3.0)));
k14=l*h;
l= x20*(rho-z20)-y20+eps1*(a21*(y10-y20)+a23*(y30-y20));
k15=l*h;
l= x20*y20-beta*z20;
k16=l*h;
l= sig*(y30-x30)+eps2*(a123*(pow(x10,2.0)*x20-pow(x30,3.0))+a123*(pow(x20,2.0)*x10-pow(x30,3.0)));
k17=l*h;
l= x30*(rho-z30)-y30+eps1*(a12*(y10-y30)+a13*(y20-y30));
k18=l*h;
l= x30*y30-beta*z30;
k19=l*h;

 
x1=x10+k11/2.0; 
y1=y10+k12/2.0;
z1=z10+k13/2.0; 

x2=x20+k14/2.0;
y2=y20+k15/2.0;
z2=z20+k16/2.0;

x3=x30+k17/2.0;
y3=y30+k18/2.0;
z3=z30+k19/2.0;
t=t0 +h/2.0;


a12=1/2 - cos(pi/3 + 2*omega*t)/3;
a13= cos(2*omega*t)/3 + 1/2;
a23=  1/2 - cos(pi/3 - 2*omega*t)/3;
a21=a12;
a31=a13;
a32=a23;
a123=1 - (2*cos(pi/3 + 2*omega*t))/3;

l= sig*(y1-x1)+eps2*(a123*(pow(x2,2.0)*x3-pow(x1,3.0))+a123*(pow(x3,2.0)*x2-pow(x1,3.0)));
k21=l*h;
l=x1*(rho-z1)-y1+eps1*(a12*(y2-y1)+a13*(y3-y1));
k22=l*h;
l=x1*y1-beta*z1;
 k23=l*h;
l= sig*(y2-x2)+eps2*(a123*(pow(x1,2.0)*x3-pow(x2,3.0))+a123*(pow(x3,2.0)*x1-pow(x2,3.0)));
k24=l*h;
l= x2*(rho-z2)-y2+eps1*(a21*(y1-y2)+a23*(y3-y2));
k25=l*h;
l= x2*y2-beta*z2;
k26=l*h;
l= sig*(y3-x3)+eps2*(a123*(pow(x1,2.0)*x2-pow(x3,3.0))+a123*(pow(x2,2.0)*x1-pow(x3,3.0)));
k27=l*h;
l= x3*(rho-z3)-y3+eps1*(a12*(y1-y3)+a13*(y2-y3));
k28=l*h;
l= x3*y3-beta*z3;
k29=l*h;
 
x1=x10+k21/2.0; 
y1=y10+k22/2.0;
z1=z10+k23/2.0; 

x2=x20+k24/2.0;
y2=y20+k25/2.0;
z2=z20+k26/2.0; 

x3=x30+k27/2.0; 
y3=y30+k28/2.0;
z3=z30+k29/2.0;
t=t0+h/2.0;


a12=1/2 - cos(pi/3 + 2*omega*t)/3;
a13= cos(2*omega*t)/3 + 1/2;
a23=  1/2 - cos(pi/3 - 2*omega*t)/3;
a21=a12;
a31=a13;
a32=a23;
a123=1 - (2*cos(pi/3 + 2*omega*t))/3;

l= sig*(y1-x1)+eps2*(a123*(pow(x2,2.0)*x3-pow(x1,3.0))+a123*(pow(x3,2.0)*x2-pow(x1,3.0)));
k31=l*h;
l=x1*(rho-z1)-y1+eps1*(a12*(y2-y1)+a13*(y3-y1));
k32=l*h;
l=x1*y1-beta*z1;
 k33=l*h;
l= sig*(y2-x2)+eps2*(a123*(pow(x1,2.0)*x3-pow(x2,3.0))+a123*(pow(x3,2.0)*x1-pow(x2,3.0)));
k34=l*h;
l= x2*(rho-z2)-y2+eps1*(a21*(y1-y2)+a23*(y3-y2));
k35=l*h;
l= x2*y2-beta*z2;
k36=l*h;
l= sig*(y3-x3)+eps2*(a123*(pow(x1,2.0)*x2-pow(x3,3.0))+a123*(pow(x2,2.0)*x1-pow(x3,3.0)));
k37=l*h;
l= x3*(rho-z3)-y3+eps1*(a12*(y1-y3)+a13*(y2-y3));
k38=l*h;
l= x3*y3-beta*z3;
k39=l*h;


x1=x10+k31; 
y1=y10+k32; 
z1=z10+k33;

x2=x20+k34;
y2=y20+k35;
z2=z20+k36; 

x3=x30+k37; 
y3=y30+k38;
z3=z30+k39;
t=t0+h;


a12=1/2 - cos(pi/3 + 2*omega*t)/3;
a13= cos(2*omega*t)/3 + 1/2;
a23=  1/2 - cos(pi/3 - 2*omega*t)/3;
a21=a12;
a31=a13;
a32=a23;
a123=1 - (2*cos(pi/3 + 2*omega*t))/3;

l= sig*(y1-x1)+eps2*(a123*(pow(x2,2.0)*x3-pow(x1,3.0))+a123*(pow(x3,2.0)*x2-pow(x1,3.0)));
k41=l*h;
l=x1*(rho-z1)-y1+eps1*(a12*(y2-y1)+a13*(y3-y1));
k42=l*h;
l=x1*y1-beta*z1;
 k43=l*h;

l= sig*(y2-x2)+eps2*(a123*(pow(x1,2.0)*x3-pow(x2,3.0))+a123*(pow(x3,2.0)*x1-pow(x2,3.0)));
k44=l*h;
l= x2*(rho-z2)-y2+eps1*(a21*(y1-y2)+a23*(y3-y2));
k45=l*h;
l= x2*y2-beta*z2;
k46=l*h;

l= sig*(y3-x3)+eps2*(a123*(pow(x1,2.0)*x2-pow(x3,3.0))+a123*(pow(x2,2.0)*x1-pow(x3,3.0)));
k47=l*h;
l= x3*(rho-z3)-y3+eps1*(a12*(y1-y3)+a13*(y2-y3));
k48=l*h;
l= x3*y3-beta*z3;
k49=l*h;
  
x10=x10+(k11+2.0*(k21+k31)+k41)/6.0;
y10=y10+(k12+2.0*(k22+k32)+k42)/6.0;
z10=z10+(k13+2.0*(k23+k33)+k43)/6.0;

x20=x20+(k14+2.0*(k24+k34)+k44)/6.0;
y20=y20+(k15+2.0*(k25+k35)+k45)/6.0;
z20=z20+(k16+2.0*(k26+k36)+k46)/6.0;

x30=x30+(k17+2.0*(k27+k37)+k47)/6.0;
y30=y30+(k18+2.0*(k28+k38)+k48)/6.0;
z30=z30+(k19+2.0*(k29+k39)+k49)/6.0;
t0=t0+h;

if(ii>tn)
{ 

 x[1]=x10;x[2]=x20;x[3]=x30;y[1]=y10;y[2]=y20;y[3]=y30;z[1]=z10;z[2]=z20;z[3]=z30;  
 fprintf(fp2,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",t0,x[1],x[2],x[3],y[1],y[2],y[3],z[1],z[2],z[3]);
  
  for (i=1;i<=3;i++)
  {
  fprintf(fp1,"%f\t%d\t%f\t%f\t%f\n",t0,i,x[i],y[i],z[i]);
  }
///fprintf(fp1,"%f\t%f\t%f\t%f\t%f\n",sig1r,sig1c,sig2r,sig2c,mle/TT);
///printf("%f\t%f\t%f\t%f\n",eps,eps2,eta,mle/TT);
}
}
}
return 0;
fclose(fp1); //// fclose(fp2);
fclose(fp2);

}

