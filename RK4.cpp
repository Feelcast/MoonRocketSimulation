#include<iostream>
#include<math.h>
#include<fstream>
#include<iomanip>
using namespace std;

//g++ -o RK4 RK4.cpp && ./RK4

double f(double t, double* xi, int n){
    double fr;
    const double om =2.6617e-6, delta  =7.014e-12, mu =0.0123;
    double r = xi[0], phi = xi[1], pr = xi[2], pp = xi[3];
    double rpr = sqrt(1+r*r-2*r*cos(phi - om*t));
    switch (n){
        case 0:
        fr = pr;
        break;
        case 1:
        fr= pp/(r*r);
        break;
        case 2:
        fr = (pr*pr)/(r*r*r)*delta*(1.0/(r*r)+(r-cos(phi - om*t))*mu/(rpr*rpr*rpr));
        break;
        case 3:
        fr=-sin(phi - om*t)*delta*mu*r/(rpr*rpr*rpr);
        break;
    }
    return fr;
}

double* xvnext(double t,double* xi,double f(double t, double* xi, int n),double h, const int n){
    double k1[n]={};
    double k2[n]={};
    double k3[n]={};
    double k4[n]={};
    double* xs = new double[n];

    for (int i = 0;i<n;i++){
        k1[i]=f(t,xi,i);
        xs[i] = xi[i] + k1[i]*h/2.0;
    }
    for (int i = 0;i<n;i++){
        k2[i] = f(t+h/2.0,xs,i);
        xs[i] = xi[i] + k2[i]*h/2.0;
    }
    for (int i = 0;i<n;i++){
        k3[i] = f(t+h/2.0,xs,i);
        xs[i] = xi[i] + k3[i]*h;
    }
    for (int i = 0;i<n;i++){
        k4[i] = f(t+h,xs,i);
    }
    for (int i = 0;i<n;i++){
        xs[i] = xi[i]+h*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
    }
    return xs;
}

void RK4sys(double ti,double* xi,double f(double t, double* xi, int n),double h, const int n, const int m){
double x[n]={};
double t = ti;
for (int i = 0;i<n;i++){
        x[i]=xi[i];
    }

ofstream file("datos.dat");
file << setprecision(8);
for (int j = 0;j<m;j++){
    file << t<<" "; 
    for (int i = 0;i<n;i++){ 
        file << x[i]<<" ";
        x[i]=xvnext(t,x,f,h,n)[i];     
    } 
    t = t +h;
    file << endl;
}
}

int main(){

const int  n =4, m=15000;
double t = 0;
double v0= 5e-5;
double r0= 0.01659;
double theta = 0.785398;
double xi[4] ={r0,-0.5,v0*cos(theta),r0*v0*sin(theta)};
double h = 2;
RK4sys(t,xi,f,h,n,m);

    return 0;
}