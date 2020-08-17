!module FFT_libraly
!  private,integer,parameter::NUM_OF_DATA=8
!  private,complex::i=(0,1.)
!  private,real::pi2=3.141592653*2.

!回転因子が複素数のFFT
program main
  implicit none
  integer,parameter::NUM_OF_DATA=8
  complex,parameter::i=(0.,1)
  real,parameter::pi2=3.141592653*2.

  integer expo_i
  integer,parameter::DATA_EXPONENT_BASE2=int(log(real(NUM_OF_DATA))/log(2.0))
  integer exponent,posi_jump, itr, dummy_itr, period_start, exp_period, count
  integer A,B,C,D

  complex f(NUM_OF_DATA),W(NUM_OF_DATA-1),temp1,temp2

  f=(0.0,0.0)
  f(1)=(1.0,0.0)


  !回転因子のテーブルを作成
  dummy_itr=1 !指定が面倒だからとりあえずイテレータを作る
  do expo_i=0,DATA_EXPONENT_BASE2-1
    do itr=0,2**(DATA_EXPONENT_BASE2-expo_i-1)-1

        W(dummy_itr) = exp(-i * pi2 * real(itr) / 2.**(DATA_EXPONENT_BASE2-expo_i) )
        dummy_itr = dummy_itr+1

    end do
  end do

  !exponentはNUM_OF_DATAの2のべき乗の指数 ex.8->2^3
  exp_period=NUM_OF_DATA/2
  period_start=1
  count=0
  do exponent=0,DATA_EXPONENT_BASE2-1
        A = ishft(1,exponent)
        B = ishft(NUM_OF_DATA,-(exponent+1))
        D = ishft(NUM_OF_DATA,-exponent)

        !posi_jumpはFFT計算の不連続点のジャンプ
        do posi_jump=0,A-1
            do itr=0,B-1

              !注意:Wのインデックスの取り方はかなり複雑
              temp1 =  f(itr+posi_jump*D+1)+f(itr+B+posi_jump*D+1)
              temp2 = ( f(itr+posi_jump*D+1)-f(itr+B+posi_jump*D+1) ) * W(period_start+count)
              f(itr +     posi_jump*D + 1)=temp1
              f(itr + B + posi_jump*D + 1)=temp2

              if(count==exp_period-1) then
                count=0
              else
                count=count+1
              end if

            end do
        end do

        period_start=period_start + exp_period
        exp_period=exp_period/2

      end do

      print *,f/NUM_OF_DATA
end program main
