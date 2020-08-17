module FFT_libraly
  implicit none
  complex,private,parameter  ::i=(0,1.)
  real,   private,parameter  ::pi2=3.141592653*2.
  integer,private            ::NUM_OF_DATA
  integer,private            ::DATA_EXPONENT_BASE2
  integer,private,allocatable::bit_reverse_index(:)
  complex,private,allocatable::W(:)
  integer,private            ::original_data_length_buf

contains
  subroutine fft_setup(nod)
    integer,intent(in)::nod
    integer expo_i,itr,dummy_itr

    original_data_length_buf = nod
    NUM_OF_DATA=2**int(log(real(nod))/log(2.0))
    DATA_EXPONENT_BASE2=int(log(real(NUM_OF_DATA))/log(2.0))
    allocate(W(NUM_OF_DATA-1))

    !回転因子のテーブルを作成
    dummy_itr=1 !指定が面倒だからとりあえずイテレータを作る
    do expo_i=0,DATA_EXPONENT_BASE2-1
      do itr=0,2**(DATA_EXPONENT_BASE2-expo_i-1)-1
          W(dummy_itr) = exp(-i * pi2 * real(itr) / 2.**(DATA_EXPONENT_BASE2-expo_i) )
          dummy_itr = dummy_itr+1
      end do
    end do

  end subroutine

  subroutine bit_reverse()
    implicit none
    integer expo,itr
    

    ! expo = NUM_OF_DATA/2
    ! bit_reverse_index(1)=0
    ! do itr=1,NUM_OF_DATA
    !   bit_reverse_index(itr)=i
    ! end do

    ! do i=1,NUM_OF_DATA
    !   do itr = i,expo
    !     bit_reverse_index(itr+expo) = NUM_OF_DATA/expo
    !   end do
    ! end do
  end subroutine

!回転因子が複素数のFFT
  subroutine FFT(output_,input_)
    complex,intent(inout)::output_(NUM_OF_DATA)
    complex,intent(in)::input_(NUM_OF_DATA)

    integer exponent,posi_jump, itr, dummy_itr, period_start, exp_period, count
    integer A,B,C,D
    complex temp1,temp2

    output_=input_

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
            temp1 =  output_(itr+posi_jump*D+1)+output_(itr+B+posi_jump*D+1)
            temp2 = (output_(itr+posi_jump*D+1)-output_(itr+B+posi_jump*D+1))*W(period_start+count)
            output_(itr+posi_jump*D+1)=temp1
            output_(itr+B+posi_jump*D+1)=temp2

            if(count==exp_period-1) then
              count=0
            else
              count=count+1
            end if

          end do
        end do

      period_start=period_start+exp_period
      exp_period=exp_period/2

    end do
  end subroutine

end module

program main
  use FFT_libraly
  implicit none
  integer,parameter::size=8
  complex out(size),in(size)
  out=(1.,0.)
  in=(0.,0.)

  call fft_setup(size)
  call bit_reverse()
  
  ! call fft(out,in)


end program
