#include<iostream>
#include<cstdlib>
#include<cmath>
#include<complex>

using namespace std;
using CD=complex<double>;
constexpr int NUM_MAX=8;
constexpr CD i(0,1.);
constexpr double pi2=3.141592653*2.;

void fft(CD F[NUM_MAX],const CD f[NUM_MAX],const int begin,const int NUM){
    if(NUM==1) return;
    CD f_[NUM_MAX];
    for(int n=begin;n<begin+NUM/2;n++)
    {
        F[n]=f[n]+f[n+NUM/2]; f_[n]=F[n];
        F[n+NUM/2]=(f[n]-f[n+NUM/2])*exp(-i*pi2*double(n-begin)/double(NUM)); f_[n+NUM/2]=F[n+NUM/2];
    }
    fft(F,f_,begin,NUM/2);
    fft(F,f_,begin+NUM/2,NUM/2);
    return;
}

int main(int argc, char* argv[])
{
    complex<double> F[NUM_MAX],f[NUM_MAX];
    f[0]=1.;
    fft(F,f,0,NUM_MAX);
    for(auto d:F)cout<<d<<" ";


    return 0;
}
