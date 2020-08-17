#include<iostream>
#include<cstdlib>
#include<cmath>
#include<complex>
#include <vector>

using namespace std;
using CD=complex<double>;
constexpr int NUM_OF_DATA=8;
constexpr CD i(0,1.);
constexpr double pi2=3.141592653*2.;

//回転因子が複素数のFFT
int main(int argc, char* argv[])
{
    CD f[NUM_OF_DATA],temp1,temp2;

    f[0]=1.; //単一パルスを入力データにする

    //回転因子のテーブルを作成
    vector<vector<CD>> W;
    W.resize(int(log2(NUM_OF_DATA)));

    for(int expo_i=0;expo_i<int(log2(NUM_OF_DATA));expo_i++){
        W[expo_i].resize(int(NUM_OF_DATA/pow(2,expo_i+1)));

        for(int itr=0;itr<W[expo_i].size();itr++){
            W[expo_i][itr]=exp(-i*pi2*double(itr*pow(2,expo_i)/NUM_OF_DATA));
        }

    }

    //exponentはNUM_OF_DATAの2のべき乗の指数 ex.8->2^3
    for(int exponent=0;exponent<int(log2(NUM_OF_DATA));exponent++){
        int A=1<<exponent;
        int B=NUM_OF_DATA>>(exponent+1);
        int D=NUM_OF_DATA>>exponent;

        //posi_jumpはFFT計算の不連続点のジャンプ
        for(int posi_jump=0;posi_jump<A;posi_jump++){
            for(int itr=0;itr<B;itr++){
                cout<<itr+posi_jump*D<<" "<<itr+B+posi_jump*D<<" "<<exponent<<endl;
                temp1= f[itr+posi_jump*D]+f[itr+B+posi_jump*D];
                temp2=(f[itr+posi_jump*D]-f[itr+B+posi_jump*D])*W[exponent][itr];
                f[itr+posi_jump*D]  =temp1;
                f[itr+B+posi_jump*D]=temp2;
            }
        }

    }
    return 0;
}
