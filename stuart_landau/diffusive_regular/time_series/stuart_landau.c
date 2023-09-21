
#include <stdio.h>  
#include <stdlib.h> 
#include <math.h> 
#include <time.h> 
int main(void)
{
srand(time(NULL));
int i,j,j1,ii,n,m,mm,j2,j3,i2,i3,i4,i5,l1,rf,rp,np,k1,IRAND,TT,k2,k3,k4,k5,nn,tn,i1,i7,na,d1,d2,N,n1;
float cum[100],v[100],gsc[100],p,temp,t0,tp0,h,l,w,sigma,rho,beta,eps,eta,RAND,ev1,alpha,a,b,c;
float ev2,znorm[100],mle,maxlyap,x10,x20,y10,y20,x30,y30,x1,x2,y1,y2,x3,y3;
float ax10,ay10,az10,ax1,ay1,az1,eps2,gc2,sig,ax20,ax2,ay20,ay2,az20,az2,ax30,ax3,ay30,ay3,az30,az3,omega;
float k14,k15,k16,k24,k25,k26,k34,k35,k36,k44,k45,k46,k42,k43,k11,k12,k13,k21,k22,k23,k31,k32,k33,k41,gc;
float tau,taus,gca,gk,gs,eca,ek,es,thetas,einh,sigr,sigc,betar,betac;
float k17,k18,k19,k27,k28,k29,k37,k38,k39,k47,k48,k49,x[4],y[4];
float sig1c,sig1r,sig2c,sig2r,ev3,alpha2,t,pi,a12,a21,a13,a31,a23,a32,a123;

FILE *fopen(),*fp1,*fp2;

fp1=fopen("omega=2.0_sig1c_sig2c_2.out","w");
fp2=fopen("L_omega=2.0_ts.out","w");
//N=100;
n=2; 
n1=2*n; //
nn=n1*(n1+1); //n*(n+1);
TT=1; 
pi=3.141592653589793;     
alpha2=2.0;
ev2=1.0;ev3=2.0; 
sigr=1.0;sigc=4.3;betar=1.0;betac=1.1;

omega=2.0;
w=sigc-betac*sigr/betar;

mm=250000; // Time Iteration
tn=245000; // Transient Time

for(i2=0;i2<=0;i2++)
{
sig1r=2.5*0.1;//2.0/100.0*i2;

for(i3=0;i3<=0;i3++)
{
sig1c=-0.5*2.5;//-5.0+10.0/100.0*i3;

for(i4=0;i4<=0;i4++)
{
sig2r=0.1*0.5;//2.0/100.0*i4;

for(i5=0;i5<=0;i5++)
{
sig2c=0.5*0.5; //-5.0+10.0/100.0*i5;



printf("%d\t%d\t%d\t%d\n",i2,i3,i4,i5);

for(i1=1;i1<=TT;i1++)
{

x10= (float)rand()/RAND_MAX; //-1.50+3.50*(float)rand()/RAND_MAX ;/// -0.07+0.05*(float)rand()/RAND_MAX; ///(float)rand()/RAND_MAX ;/// (-1.5+rand()*3.50)*0.01
y10= (float)rand()/RAND_MAX; //-7.0+8.0*(float)rand()/RAND_MAX ;///0.12*(float)rand()/RAND_MAX; (-7.0+rand()*8.0)*0.01
x20= (float)rand()/RAND_MAX; //-1.50+3.50*(float)rand()/RAND_MAX ;/// -0.07+0.05*(float)rand()/RAND_MAX; ///(float)rand()/RAND_MAX ;/// (-1.5+rand()*3.50)*0.01
y20= (float)rand()/RAND_MAX; //-7.0+8.0*(float)rand()/RAND_MAX ;///0.12*(float)rand()/RAND_MAX; (-7.0+rand()*8.0)*0.01
x30= (float)rand()/RAND_MAX; //-1.50+3.50*(float)rand()/RAND_MAX ;/// -0.07+0.05*(float)rand()/RAND_MAX; ///(float)rand()/RAND_MAX ;/// (-1.5+rand()*3.50)*0.01
y30= (float)rand()/RAND_MAX; //-7.0+8.0*(float)rand()/RAND_MAX ;///0.12*(float)rand()/RAND_MAX; (-7.0+rand()*8.0)*0.01




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

 l= sigr*x10-sigc*y10-betar*x10*(pow(x10,2.0)+pow(y10,2.0))+betac*y10*(pow(x10,2.0)+pow(y10,2.0))+sig1r*a12*(x20-x10)-sig1c*a12*(y20-y10)+sig1r*a13*(x30-x10)-sig1c*a13*(y30-y10) \
    +2*sig2r*a123*(x20*x30-y20*y30+pow(x10,2.0)-pow(y10,2.0))-2*sig2c*a123*(x20*y30+x30*y20-2*x10*y10);
k11=l*h;
l=sigr*y10+sigc*x10-betac*x10*(pow(x10,2.0)+pow(y10,2.0))-betar*y10*(pow(x10,2.0)+pow(y10,2.0))+sig1c*a12*(x20-x10)+sig1r*a12*(y20-y10)+sig1c*a13*(x30-x10)+sig1r*a13*(y30-y10) \
    +2*sig2c*a123*(x20*x30-y20*y30+pow(x10,2.0)-pow(y10,2.0))+2*sig2r*a123*(x20*y30+x30*y20-2*x10*y10);
k12=l*h;

 l= sigr*x20-sigc*y20-betar*x20*(pow(x20,2.0)+pow(y20,2.0))+betac*y20*(pow(x20,2.0)+pow(y20,2.0))+sig1r*a21*(x10-x20)-sig1c*a21*(y10-y20)+sig1r*a23*(x30-x20)-sig1c*a23*(y30-y20) \
    +2*sig2r*a123*(x10*x30-y10*y30+pow(x20,2.0)-pow(y20,2.0))-2*sig2c*a123*(x10*y30+x10*y30-2*x20*y20);
k13=l*h;
l=sigr*y20+sigc*x20-betac*x20*(pow(x20,2.0)+pow(y20,2.0))-betar*y20*(pow(x20,2.0)+pow(y20,2.0))+sig1c*a21*(x10-x20)+sig1r*a21*(y10-y20)+sig1c*a23*(x30-x20)+sig1r*a23*(y30-y20) \
    +2*sig2c*a123*(x10*x30-y10*y30+pow(x20,2.0)-pow(y20,2.0))+2*sig2r*a123*(x10*y30+x10*y30-2*x20*y20);
k14=l*h;

 l= sigr*x30-sigc*y30-betar*x30*(pow(x30,2.0)+pow(y30,2.0))+betac*y30*(pow(x30,2.0)+pow(y30,2.0))+sig1r*a31*(x10-x30)-sig1c*a31*(y10-y30)+sig1r*a32*(x20-x30)-sig1c*a32*(y20-y30) \
    +2*sig2r*a123*(x10*x20-y10*y20+pow(x30,2.0)-pow(y30,2.0))-2*sig2c*a123*(x10*y20+x10*y20-2*x30*y30);
k15=l*h;
l=sigr*y30+sigc*x30-betac*x30*(pow(x30,2.0)+pow(y30,2.0))-betar*y30*(pow(x30,2.0)+pow(y30,2.0))+sig1c*a31*(x10-x30)+sig1r*a31*(y10-y30)+sig1c*a32*(x20-x30)+sig1r*a32*(y20-y30) \
    +2*sig2c*a123*(x10*x20-y10*y20+pow(x30,2.0)-pow(y30,2.0))+2*sig2r*a123*(x10*y20+x10*y20-2*x30*y30);
k16=l*h;

 
x1=x10+k11/2.0; 
y1=y10+k12/2.0; 
x2=x20+k13/2.0;
y2=y20+k14/2.0;
x3=x30+k15/2.0;
y3=y30+k16/2.0;
t=t0 +h/2.0;


a12=1/2 - cos(pi/3 + 2*omega*t)/3;
a13= cos(2*omega*t)/3 + 1/2;
a23=  1/2 - cos(pi/3 - 2*omega*t)/3;
a21=a12;
a31=a13;
a32=a23;
a123=1 - (2*cos(pi/3 + 2*omega*t))/3;

 l= sigr*x1-sigc*y1-betar*x1*(pow(x1,2.0)+pow(y1,2.0))+betac*y1*(pow(x1,2.0)+pow(y1,2.0))+sig1r*a12*(x2-x1)-sig1c*a12*(y2-y1)+sig1r*a13*(x3-x1)-sig1c*a13*(y3-y1) \
    +2*sig2r*a123*(x2*x3-y2*y3+pow(x1,2.0)-pow(y1,2.0))-2*sig2c*a123*(x2*y3+x3*y2-2*x1*y1);
k21=l*h;
l=sigr*y1+sigc*x1-betac*x1*(pow(x1,2.0)+pow(y1,2.0))-betar*y1*(pow(x1,2.0)+pow(y1,2.0))+sig1c*a12*(x2-x1)+sig1r*a12*(y2-y1)+sig1c*a13*(x3-x1)+sig1r*a13*(y3-y1) \
    +2*sig2c*a123*(x2*x3-y2*y3+pow(x1,2.0)-pow(y1,2.0))+2*sig2r*a123*(x2*y3+x3*y2-2*x1*y1);
k22=l*h;

 l= sigr*x2-sigc*y2-betar*x2*(pow(x2,2.0)+pow(y2,2.0))+betac*y2*(pow(x2,2.0)+pow(y2,2.0))+sig1r*a21*(x1-x2)-sig1c*a21*(y1-y2)+sig1r*a23*(x3-x2)-sig1c*a23*(y3-y2) \
    +2*sig2r*a123*(x1*x3-y1*y3+pow(x2,2.0)-pow(y2,2.0))-2*sig2c*a123*(x1*y3+x1*y3-2*x2*y2);
k23=l*h;
l=sigr*y2+sigc*x2-betac*x2*(pow(x2,2.0)+pow(y2,2.0))-betar*y2*(pow(x2,2.0)+pow(y2,2.0))+sig1c*a21*(x1-x2)+sig1r*a21*(y1-y2)+sig1c*a23*(x3-x2)+sig1r*a23*(y3-y2) \
    +2*sig2c*a123*(x1*x3-y1*y3+pow(x2,2.0)-pow(y2,2.0))+2*sig2r*a123*(x1*y3+x1*y3-2*x2*y2);
k24=l*h;

 l= sigr*x3-sigc*y3-betar*x3*(pow(x3,2.0)+pow(y3,2.0))+betac*y3*(pow(x3,2.0)+pow(y3,2.0))+sig1r*a31*(x1-x3)-sig1c*a31*(y1-y3)+sig1r*a32*(x2-x3)-sig1c*a32*(y2-y3) \
    +2*sig2r*a123*(x1*x2-y1*y2+pow(x3,2.0)-pow(y3,2.0))-2*sig2c*a123*(x1*y2+x1*y2-2*x3*y3);
k25=l*h;
l=sigr*y3+sigc*x3-betac*x3*(pow(x3,2.0)+pow(y3,2.0))-betar*y3*(pow(x3,2.0)+pow(y3,2.0))+sig1c*a31*(x1-x3)+sig1r*a31*(y1-y3)+sig1c*a32*(x2-x3)+sig1r*a32*(y2-y3) \
    +2*sig2c*a123*(x1*x2-y1*y2+pow(x3,2.0)-pow(y3,2.0))+2*sig2r*a123*(x1*y2+x1*y2-2*x3*y3);
k26=l*h;
 
x1=x10+k21/2.0; 
y1=y10+k22/2.0; 
x2=x20+k23/2.0;
y2=y20+k24/2.0; 
x3=x30+k25/2.0; 
y3=y30+k26/2.0;
t=t0+h/2.0;


a12=1/2 - cos(pi/3 + 2*omega*t)/3;
a13= cos(2*omega*t)/3 + 1/2;
a23=  1/2 - cos(pi/3 - 2*omega*t)/3;
a21=a12;
a31=a13;
a32=a23;
a123=1 - (2*cos(pi/3 + 2*omega*t))/3;

 l= sigr*x1-sigc*y1-betar*x1*(pow(x1,2.0)+pow(y1,2.0))+betac*y1*(pow(x1,2.0)+pow(y1,2.0))+sig1r*a12*(x2-x1)-sig1c*a12*(y2-y1)+sig1r*a13*(x3-x1)-sig1c*a13*(y3-y1) \
    +2*sig2r*a123*(x2*x3-y2*y3+pow(x1,2.0)-pow(y1,2.0))-2*sig2c*a123*(x2*y3+x3*y2-2*x1*y1);
k31=l*h;
l=sigr*y1+sigc*x1-betac*x1*(pow(x1,2.0)+pow(y1,2.0))-betar*y1*(pow(x1,2.0)+pow(y1,2.0))+sig1c*a12*(x2-x1)+sig1r*a12*(y2-y1)+sig1c*a13*(x3-x1)+sig1r*a13*(y3-y1) \
    +2*sig2c*a123*(x2*x3-y2*y3+pow(x1,2.0)-pow(y1,2.0))+2*sig2r*a123*(x2*y3+x3*y2-2*x1*y1);
k32=l*h;

 l= sigr*x2-sigc*y2-betar*x2*(pow(x2,2.0)+pow(y2,2.0))+betac*y2*(pow(x2,2.0)+pow(y2,2.0))+sig1r*a21*(x1-x2)-sig1c*a21*(y1-y2)+sig1r*a23*(x3-x2)-sig1c*a23*(y3-y2) \
    +2*sig2r*a123*(x1*x3-y1*y3+pow(x2,2.0)-pow(y2,2.0))-2*sig2c*a123*(x1*y3+x1*y3-2*x2*y2);
k33=l*h;
l=sigr*y2+sigc*x2-betac*x2*(pow(x2,2.0)+pow(y2,2.0))-betar*y2*(pow(x2,2.0)+pow(y2,2.0))+sig1c*a21*(x1-x2)+sig1r*a21*(y1-y2)+sig1c*a23*(x3-x2)+sig1r*a23*(y3-y2) \
    +2*sig2c*a123*(x1*x3-y1*y3+pow(x2,2.0)-pow(y2,2.0))+2*sig2r*a123*(x1*y3+x1*y3-2*x2*y2);
k34=l*h;

 l= sigr*x3-sigc*y3-betar*x3*(pow(x3,2.0)+pow(y3,2.0))+betac*y3*(pow(x3,2.0)+pow(y3,2.0))+sig1r*a31*(x1-x3)-sig1c*a31*(y1-y3)+sig1r*a32*(x2-x3)-sig1c*a32*(y2-y3) \
    +2*sig2r*a123*(x1*x2-y1*y2+pow(x3,2.0)-pow(y3,2.0))-2*sig2c*a123*(x1*y2+x1*y2-2*x3*y3);
k35=l*h;
l=sigr*y3+sigc*x3-betac*x3*(pow(x3,2.0)+pow(y3,2.0))-betar*y3*(pow(x3,2.0)+pow(y3,2.0))+sig1c*a31*(x1-x3)+sig1r*a31*(y1-y3)+sig1c*a32*(x2-x3)+sig1r*a32*(y2-y3) \
    +2*sig2c*a123*(x1*x2-y1*y2+pow(x3,2.0)-pow(y3,2.0))+2*sig2r*a123*(x1*y2+x1*y2-2*x3*y3);
k36=l*h;

x1=x10+k31; 
y1=y10+k32; 
x2=x20+k33;
y2=y20+k34; 
x3=x30+k35; 
y3=y30+k36;
t=t0+h;


a12=1/2 - cos(pi/3 + 2*omega*t)/3;
a13= cos(2*omega*t)/3 + 1/2;
a23=  1/2 - cos(pi/3 - 2*omega*t)/3;
a21=a12;
a31=a13;
a32=a23;
a123=1 - (2*cos(pi/3 + 2*omega*t))/3;

 l= sigr*x1-sigc*y1-betar*x1*(pow(x1,2.0)+pow(y1,2.0))+betac*y1*(pow(x1,2.0)+pow(y1,2.0))+sig1r*a12*(x2-x1)-sig1c*a12*(y2-y1)+sig1r*a13*(x3-x1)-sig1c*a13*(y3-y1) \
    +2*sig2r*a123*(x2*x3-y2*y3+pow(x1,2.0)-pow(y1,2.0))-2*sig2c*a123*(x2*y3+x3*y2-2*x1*y1);
k41=l*h;
l=sigr*y1+sigc*x1-betac*x1*(pow(x1,2.0)+pow(y1,2.0))-betar*y1*(pow(x1,2.0)+pow(y1,2.0))+sig1c*a12*(x2-x1)+sig1r*a12*(y2-y1)+sig1c*a13*(x3-x1)+sig1r*a13*(y3-y1) \
    +2*sig2c*a123*(x2*x3-y2*y3+pow(x1,2.0)-pow(y1,2.0))+2*sig2r*a123*(x2*y3+x3*y2-2*x1*y1);
k42=l*h;

 l= sigr*x2-sigc*y2-betar*x2*(pow(x2,2.0)+pow(y2,2.0))+betac*y2*(pow(x2,2.0)+pow(y2,2.0))+sig1r*a21*(x1-x2)-sig1c*a21*(y1-y2)+sig1r*a23*(x3-x2)-sig1c*a23*(y3-y2) \
    +2*sig2r*a123*(x1*x3-y1*y3+pow(x2,2.0)-pow(y2,2.0))-2*sig2c*a123*(x1*y3+x1*y3-2*x2*y2);
k43=l*h;
l=sigr*y2+sigc*x2-betac*x2*(pow(x2,2.0)+pow(y2,2.0))-betar*y2*(pow(x2,2.0)+pow(y2,2.0))+sig1c*a21*(x1-x2)+sig1r*a21*(y1-y2)+sig1c*a23*(x3-x2)+sig1r*a23*(y3-y2) \
    +2*sig2c*a123*(x1*x3-y1*y3+pow(x2,2.0)-pow(y2,2.0))+2*sig2r*a123*(x1*y3+x1*y3-2*x2*y2);
k44=l*h;

 l= sigr*x3-sigc*y3-betar*x3*(pow(x3,2.0)+pow(y3,2.0))+betac*y3*(pow(x3,2.0)+pow(y3,2.0))+sig1r*a31*(x1-x3)-sig1c*a31*(y1-y3)+sig1r*a32*(x2-x3)-sig1c*a32*(y2-y3) \
    +2*sig2r*a123*(x1*x2-y1*y2+pow(x3,2.0)-pow(y3,2.0))-2*sig2c*a123*(x1*y2+x1*y2-2*x3*y3);
k45=l*h;
l=sigr*y3+sigc*x3-betac*x3*(pow(x3,2.0)+pow(y3,2.0))-betar*y3*(pow(x3,2.0)+pow(y3,2.0))+sig1c*a31*(x1-x3)+sig1r*a31*(y1-y3)+sig1c*a32*(x2-x3)+sig1r*a32*(y2-y3) \
    +2*sig2c*a123*(x1*x2-y1*y2+pow(x3,2.0)-pow(y3,2.0))+2*sig2r*a123*(x1*y2+x1*y2-2*x3*y3);
k46=l*h;

  
x10=x10+(k11+2.0*(k21+k31)+k41)/6.0;
y10=y10+(k12+2.0*(k22+k32)+k42)/6.0;
x20=x20+(k13+2.0*(k23+k33)+k43)/6.0;
y20=y20+(k14+2.0*(k24+k34)+k44)/6.0;
x30=x30+(k15+2.0*(k25+k35)+k45)/6.0;
y30=y30+(k16+2.0*(k26+k36)+k46)/6.0;
t0=t0+h;

if(ii>tn)
{ 

 x[1]=x10;x[2]=x20;x[3]=x30;y[1]=y10;y[2]=y20;y[3]=y30;  
 fprintf(fp2,"%f\t%f\t%f\t%f\t%f\t%f\t%f\n",t0,x[1],x[2],x[3],y[1],y[2],y[3]);
  
  for (i=1;i<=3;i++)
  {
  fprintf(fp1,"%f\t%d\t%f\t%f\n",t0,i,x[i],y[i]);
  }
///fprintf(fp1,"%f\t%f\t%f\t%f\t%f\n",sig1r,sig1c,sig2r,sig2c,mle/TT);
///printf("%f\t%f\t%f\t%f\n",eps,eps2,eta,mle/TT);
}
}
}
}
}
}
}
return 0;
fclose(fp1); //// fclose(fp2);
fclose(fp2);

}

