#include <stdio.h>  
#include <stdlib.h> 
#include <math.h> 
#include <time.h> 
int main(void)
{
 srand(1);  // srand(time(NULL)); //
int i,j,j1,ii,n,m,mm,j2,j3,i2,i3,i4,i5,l1,rf,rp,np,k1,IRAND,TT,k2,k3,k4,k5,nn,tn,i1,i7,na,d1,d2,N,n1;
float cum[300000],v[300000],gsc[300000],p,temp,t0,tp0,h,l,w,sigma,rho,beta,eps,eta,RAND,alpha;
float znorm[300000],mle,maxlyap,x10,v20,y10,n20,z10,s20,x1,v2,y1,n2,z1,s2;
float eps1,eps2,sig,omega,sigr,sigc,betar,betac;
//float k14,k15,k16,k24,k25,k26,k34,k35,k36,k44,k45,k46,k42,k43,k11,k12,k13,k21,k22,k23,k31,k32,k33,k41,gc;
float k[4][301],sum1,sum2,ax0[101],ax[101],ay0[101],ay[101],random_value;
float sig1c,sig1r,sig2c,sig2r,ev3,alpha2,t,C[101][101],ev[101];

FILE *fopen(),*fp1; ///,*fp2;

fp1=fopen("N=50_eps1_eps2=1.0_5.out","w");
///fp2=fopen("L_ts.out","w");
//N=100;
n=(100-1); 
n1=2*n; //
nn=n1*(n1+1); //n*(n+1);
TT=1; 
alpha2=2.0;
 
sigr=1.0;sigc=4.3;betar=1.0;betac=1.1;
sig1r=0.1;sig2r=0.1;sig1c=-0.5;sig2c=0.5;
//omega=0.8;
w=sigc-betac*sigr/betar;

mm=200000;//250000; // Time Iteration
tn=0;//200000; // Transient Time


///////////////////////////

for(j=1;j<=n;j++)
{
 ev[j]=10.0*(float)rand()/RAND_MAX;
}

for(i=1;i<=n;i++)
{
for(j=i+1;j<=n;j++)
{
random_value= (float)rand()/RAND_MAX;
C[i][j]= random_value;
C[j][i]=-random_value;
}
}

for(i=1;i<=n;i++)
{
C[i][i]=0.0;
}
////////////////////////////////


for(i2=0;i2<=50;i2++)
{
eps1=10.0/50.0*i2;

for(i3=0;i3<=0;i3++)
{
eps2=1.0; //5.0/30.0*i3;


////////////////////////////
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
for(j=1;j<=n;j++)
{
ax0[j]=v[(2*j-1)*n1+1+j1]; 
ay0[j]=v[(2*j)*n1+1+j1];
} 
for(j=1;j<=n;j++)
{
sum1=0.0;
sum2=0.0;
for(m=1;m<=n;m++)
{
sum1=sum1+C[j][m]*ax0[m];
sum2=sum2+C[j][m]*ay0[m];
}
l= sum1-2*sigr*ax0[j]-ev[j]*(eps1*sig1r*ax0[j]-eps1*sig1c*ay0[j]+2*alpha2*sqrt(sigr/betar)*((cos(w*t0)*eps2*sig2r-sin(w*t0)*eps2*sig2c)*ax0[j]+(-cos(w*t0)*eps2*sig2c-sin(w*t0)*eps2*sig2r)*ay0[j])); 
k[1][2*j-1]=l*h;
l=sum2-2*betac*(sigr/betar)*ax0[j]-ev[j]*(eps1*sig1c*ax0[j]+eps1*sig1r*ay0[j]+2*alpha2*sqrt(sigr/betar)*((sin(w*t0)*eps2*sig2r+cos(w*t0)*eps2*sig2c)*ax0[j]\
+(-sin(w*t0)*eps2*sig2c+cos(w*t0)*eps2*sig2r)*ay0[j])); 
 k[1][2*j]=l*h;
}

for(j=1;j<=n;j++)
{
ax[j]=ax0[j]+k[1][2*j-1]/2.0; 
ay[j]=ay0[j]+k[1][2*j]/2.0;
} 
t=t0 +h/2.0;

for(j=1;j<=n;j++)
{
sum1=0.0;
sum2=0.0;
for(m=1;m<=n;m++)
{
sum1=sum1+C[j][m]*ax[m];
sum2=sum2+C[j][m]*ay[m];
}
l= sum1-2*sigr*ax[j]-ev[j]*(eps1*sig1r*ax[j]-eps1*sig1c*ay[j]+2*alpha2*sqrt(sigr/betar)*((cos(w*t)*eps2*sig2r-sin(w*t)*eps2*sig2c)*ax[j]+(-cos(w*t)*eps2*sig2c-sin(w*t)*eps2*sig2r)*ay[j])); 
k[2][2*j-1]=l*h;
l=sum2-2*betac*(sigr/betar)*ax[j]-ev[j]*(eps1*sig1c*ax[j]+eps1*sig1r*ay[j]+2*alpha2*sqrt(sigr/betar)*((sin(w*t)*eps2*sig2r+cos(w*t)*eps2*sig2c)*ax[j]\
+(-sin(w*t)*eps2*sig2c+cos(w*t)*eps2*sig2r)*ay[j])); 
 k[2][2*j]=l*h;
}

for(j=1;j<=n;j++)
{
ax[j]=ax0[j]+k[2][2*j-1]/2.0;
ay[j]=ay0[j]+k[2][2*j]/2.0;
}
t=t0+h/2.0;

for(j=1;j<=n;j++)
{
sum1=0.0;
sum2=0.0;
for(m=1;m<=n;m++)
{
sum1=sum1+C[j][m]*ax[m];
sum2=sum2+C[j][m]*ay[m];
}
l= sum1-2*sigr*ax[j]-ev[j]*(eps1*sig1r*ax[j]-eps1*sig1c*ay[j]+2*alpha2*sqrt(sigr/betar)*((cos(w*t)*eps2*sig2r-sin(w*t)*eps2*sig2c)*ax[j]+(-cos(w*t)*eps2*sig2c-sin(w*t)*eps2*sig2r)*ay[j])); 
k[3][2*j-1]=l*h;
l=sum2-2*betac*(sigr/betar)*ax[j]-ev[j]*(eps1*sig1c*ax[j]+eps1*sig1r*ay[j]+2*alpha2*sqrt(sigr/betar)*((sin(w*t)*eps2*sig2r+cos(w*t)*eps2*sig2c)*ax[j]\
+(-sin(w*t)*eps2*sig2c+cos(w*t)*eps2*sig2r)*ay[j])); 
 k[3][2*j]=l*h;
}

for(j=1;j<=n;j++)
{
ax[j]=ax0[j]+k[3][2*j-1];
ay[j]=ay0[j]+k[3][2*j];
}
t=t0+h;

for(j=1;j<=n;j++)
{
sum1=0.0;
sum2=0.0;
for(m=1;m<=n;m++)
{
sum1=sum1+C[j][m]*ax[m];
sum2=sum2+C[j][m]*ay[m];
}
l= sum1-2*sigr*ax[j]-ev[j]*(eps1*sig1r*ax[j]-eps1*sig1c*ay[j]+2*alpha2*sqrt(sigr/betar)*((cos(w*t)*eps2*sig2r-sin(w*t)*eps2*sig2c)*ax[j]+(-cos(w*t)*eps2*sig2c-sin(w*t)*eps2*sig2r)*ay[j])); 
k[4][2*j-1]=l*h;
l=sum2-2*betac*(sigr/betar)*ax[j]-ev[j]*(eps1*sig1c*ax[j]+eps1*sig1r*ay[j]+2*alpha2*sqrt(sigr/betar)*((sin(w*t)*eps2*sig2r+cos(w*t)*eps2*sig2c)*ax[j]\
+(-sin(w*t)*eps2*sig2c+cos(w*t)*eps2*sig2r)*ay[j])); 
 k[4][2*j]=l*h;
}

for(j=1;j<=n;j++)
{
ax0[j]=ax0[j]+(k[1][2*j-1]+2.0*(k[2][2*j-1]+k[3][2*j-1])+k[4][2*j-1])/6.0;
ay0[j]=ay0[j]+(k[1][2*j]+2.0*(k[2][2*j]+k[3][2*j])+k[4][2*j])/6.0;
}
t0=t0+h;

for(j=1;j<=n;j++)
{
v[(2*j-1)*n1+1+j1]=ax0[j]; 
v[(2*j)*n1+1+j1]=ay0[j]; 
}

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
fprintf(fp1,"%f\t%f\t%f\n",eps1,eps2,mle/TT);
///printf("%f\t%f\t%f\t%f\n",eps,eps2,eta,mle/TT);
}
}

return 0;
fclose(fp1); //// fclose(fp2);

}


