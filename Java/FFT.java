import java.lang.Math;

//回転因子が複素数のFFT
public class FFT
{
    private int NUM_OF_DATA;
    private final double pi2=3.141592653*2.;
    
    private double[] W_sin,W_cos;

    private int DATA_EXPONENT_BASE2;

    public FFT(int NUM_OF_DATA_){
        NUM_OF_DATA=NUM_OF_DATA_;

        //回転因子のテーブルを作成
        W_sin=new double[NUM_OF_DATA-1];
        W_cos=new double[NUM_OF_DATA-1];

        DATA_EXPONENT_BASE2=(int)( Math.log((double)NUM_OF_DATA)/Math.log(2.0) );
    
    int dummy_itr=0; //指定が面倒だからとりあえずイテレータを作る
    for(int expo_i=0; expo_i<DATA_EXPONENT_BASE2; expo_i++){
        for(int itr=0; itr<Math.pow(2,DATA_EXPONENT_BASE2-expo_i-1); itr++){
        W_cos[dummy_itr] = Math.cos( -pi2 * (double)itr / Math.pow(2.,DATA_EXPONENT_BASE2-expo_i) );
        W_sin[dummy_itr] = Math.sin( -pi2 * (double)itr / Math.pow(2.,DATA_EXPONENT_BASE2-expo_i) );
        dummy_itr++;

      }
    }

    }

    public double power_spector(double[] data_real){
        double[] data_imag=new double[data_real.length];

        int exp_period=NUM_OF_DATA/2, period_start=0, count=0;

        //exponentはNUM_OF_DATAの2のべき乗の指数 ex.8->2^3
        for(int exponent=0;exponent<DATA_EXPONENT_BASE2;exponent++){
            int A=1<<exponent;
            int B=NUM_OF_DATA>>(exponent+1);
            int D=NUM_OF_DATA>>exponent;

            double temp1_r,temp1_i,temp2_r,temp2_i;
            //posi_jumpはFFT計算の不連続点のジャンプ
            for(int posi_jump=0;posi_jump<A;posi_jump++){
                for(int itr=0;itr<B;itr++){                    
                    temp1_r = data_real[itr+posi_jump*D]+data_real[itr+B+posi_jump*D];
                    temp1_i = data_imag[itr+posi_jump*D]+data_imag[itr+B+posi_jump*D];

                    temp2_r = data_real[itr+posi_jump*D]-data_real[itr+B+posi_jump*D];
                    temp2_i = data_imag[itr+posi_jump*D]-data_imag[itr+B+posi_jump*D];

                    data_real[itr+posi_jump*D]  = temp1_r;
                    data_imag[itr+posi_jump*D]  = temp1_i;
                    
                    data_real[itr+B+posi_jump*D] = temp2_r*W_cos[period_start+count] + temp2_i*W_sin[period_start+count];
                    data_imag[itr+B+posi_jump*D] = temp2_r*W_sin[period_start+count] + temp2_i*W_cos[period_start+count];

                    count=(count==exp_period-1) ? 0 : count+1;

                }
            }
            period_start=period_start+exp_period;
            exp_period=exp_period/2;
        }
        double power=0.0;

        //ユニタリのFFTにするため各要素にsqrt(NUM_OF_DATA)で割る。->2乗総和の最後にNUM_OF_DATAで割る。
        //ユニタリにするとパワースペクトルの総和が元の関数のパワースペクトルの総和に一致する
        for(int d_i=0;d_i<data_real.length;d_i++){
            power += Math.pow(data_real[d_i],2) + Math.pow(data_imag[d_i],2);
        }
        return power/NUM_OF_DATA;
    }
}

