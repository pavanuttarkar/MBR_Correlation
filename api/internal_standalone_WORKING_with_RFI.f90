subroutine external_loop (comf_X, comf_Y, comf1_X, comf1_Y, avg, le, fftsize, RFI_X1, RFI_Y1, RFI_X2, RFI_Y2, &
                              &  len_1_x, len_1_y, len_2_x, len_2_y, creal1, creal2, creal3, creal4,  &
                              &  creal8, creal9)
  
  !f2py threadsafe
  integer(kind=4), intent(in):: len_1_x
  integer(kind=4), intent(in):: len_2_x
  integer(kind=4), intent(in):: len_1_y
  integer(kind=4), intent(in):: len_2_y

  integer(kind=4), intent(in):: fftsize
  integer(kind=4), intent(in):: avg
  integer(kind=4), intent(in):: le
  
  integer(kind=1), intent(in), dimension(len_1_x)::comf_X
  integer(kind=1), intent(in), dimension(len_1_y)::comf_Y
  integer(kind=1), intent(in), dimension(len_2_x)::comf1_X
  integer(kind=1), intent(in), dimension(len_2_y)::comf1_Y

  !RFI Matrix!
  integer(kind=1), intent(in), dimension(fftsize)::RFI_X1
  integer(kind=1), intent(in), dimension(fftsize)::RFI_Y1
  integer(kind=1), intent(in), dimension(fftsize)::RFI_X2
  integer(kind=1), intent(in), dimension(fftsize)::RFI_Y2


  !Auto Correlation!
  real             , intent(out), dimension(fftsize, le)::creal1
  real             , intent(out), dimension(fftsize, le)::creal2
  real             , intent(out), dimension(fftsize, le)::creal3
  real             , intent(out), dimension(fftsize, le)::creal4

  !Cross Pol Correlation!
  !double complex, intent(out), dimension(fftsize, le)::creal5
  !double complex, intent(out), dimension(fftsize, le)::creal6

  !Cross Pol leakage!
  !double complex, intent(out), dimension(fftsize, le)::pollkg1
  !double complex, intent(out), dimension(fftsize, le)::pollkg2

  !Cross Correlation!
  double complex, intent(out), dimension(fftsize, le)::creal8
  double complex, intent(out), dimension(fftsize, le)::creal9
 
  integer(kind=4) :: p_c

  integer(kind=8)  ::fir
  integer(kind=8)  ::sed
  integer(kind=8)  ::sed1
  integer(kind=8)  ::fir1
  integer(kind=8)  ::m

  integer(kind=1), dimension(fftsize*avg*2)::chunk_X
  integer(kind=1), dimension(fftsize*avg*2):: chunk_Y
  integer(kind=1), dimension(fftsize*avg*2):: chunk1_X
  integer(kind=1), dimension(fftsize*avg*2):: chunk1_Y

  real   (kind=8), dimension(fftsize*2)  ::X1 
  double complex, dimension(fftsize*2)  ::fX1


  double complex, dimension(fftsize)  ::fX1_1
  double complex, dimension(fftsize)  ::fY1_1
  double complex, dimension(fftsize)  ::fX2_1
  double complex, dimension(fftsize)  ::fY2_1
  
  real :: start, finish, start1, finish1, start2, finish2, cpu_t
  real(kind=8)  , dimension(fftsize)  ::ONE  
  real(kind=8)  , dimension(fftsize)  ::TWO  
  real(kind=8)  , dimension(fftsize)  ::THREE
  real(kind=8)  , dimension(fftsize)  ::FOUR 
  !double complex, dimension(fftsize)  ::FIVE 
  !double complex, dimension(fftsize)  ::SIX  
  double complex, dimension(fftsize)  ::SEVEN
  double complex, dimension(fftsize)  ::EIGHT
  !double complex, dimension(fftsize)  ::NINE 
  !double complex, dimension(fftsize)  ::TEN  

  integer(kind=8)::plan_forward1!, plan_forward2, plan_forward3, plan_forward4
  real(kind=8)::sqt, fftsize1!iret
  !call dfftw_init_threads(iret)
  !print *, iret
  !call dfftw_plan_with_nthreads(8)
  !sed1  =  avg*fftsize*2
  avg1    = avg-1
  cpu_t=0
  p_c=1
  m=1
  le1 = le-1
  fftsize1=fftsize
  sqt=sqrt(fftsize1)
  print *, 'Planning FFT'
  call dfftw_plan_dft_r2c_1d ( plan_forward1, 512, X1, fX1, FFTW_EXHAUSTIVE )!call dfftw_plan_dft_r2c_1d_ ( plan_forward1, 512, X1, fX1, FFTW_ESTIMATE )
  ONE   =       0.0!(/ (0, I = 1,fftsize) /)!0
  TWO   =       0.0!(/ (0, I = 1,fftsize) /)!0
  THREE =       0.0!(/ (0, I = 1,fftsize) /)!0
  FOUR  =       0.0!(/ (0, I = 1,fftsize) /)!0
  !FIVE  =       (/ (0, I = 1,fftsize) /)!(0,0)
  !SIX   =       (/ (0, I = 1,fftsize) /)!(0, 0)
  SEVEN =       (0.0, 0.0)!(/ (0, I = 1,fftsize) /)!(0,0)
  EIGHT =       (0.0, 0.0)!(/ (0, I = 1,fftsize) /)!(0,0)
  !NINE  =       (/ (0, I = 1,fftsize) /)!(0,0)
  !TEN   =       (/ (0, I = 1,fftsize) /)!(0,0)
  !print *, ONE
  !print *, le
  call cpu_time(start)
  do p_c =1, le 
    !print *, p_c
    fir             = avg*(p_c-1)*fftsize*2+1
    sed             = avg*(p_c)*fftsize*2

    chunk_X         = comf_X(fir:sed)
    chunk1_X        = comf1_X(fir:sed)
    chunk_Y         = comf_Y(fir:sed)
    chunk1_Y        = comf1_Y(fir:sed)
    !call cpu_time(start1)
    m=1
    do m  =1, avg
      fir1=  (m-1)*fftsize*2+1
      sed1=  (m)*fftsize*2
      
      
      X1  =  chunk_X(fir1:sed1)
      call dfftw_execute( plan_forward1 )!call dfftw_execute_ ( plan_forward1 )
      fX1_1 = fX1(1:256)*RFI_X1/sqt
      ONE  = ONE + fX1_1*CONJG(fX1_1)
      !print *, fX1_1!ONE
      
      X1  =  chunk_Y(fir1:sed1)
      call dfftw_execute( plan_forward1 ) !call dfftw_execute_ ( plan_forward2 )
      fY1_1 = fX1(1:256)*RFI_Y1/sqt
      TWO      = TWO + fY1_1*CONJG(fY1_1)

      X1  =  chunk1_X(fir1:sed1)
      call dfftw_execute( plan_forward1 ) !call dfftw_execute_ ( plan_forward3 )
      fX2_1 = fX1(1:256)*RFI_X2/sqt
      THREE=THREE + fX2_1*CONJG(fX2_1)

      X1  =  chunk1_Y(fir1:sed1)
      call dfftw_execute( plan_forward1 ) !call dfftw_execute_ ( plan_forward4 )
      fY2_1 = fX1(1:256)*RFI_Y2/sqt
      FOUR=FOUR  + fY2_1*CONJG(fY2_1)


      !FIVE         = FIVE  + CONJG(fX1_1)*fY1_1
      !SIX          = SIX   + CONJG(fX2_1)*fY2_1
      SEVEN        = SEVEN + fX1_1*CONJG(fX2_1)
      EIGHT        = EIGHT + fY1_1*CONJG(fY2_1)
      !NINE         = NINE  + fY1_1*CONJG(fY2_1)
      !TEN          = TEN   + CONJG(fX2_1)*fY1_1

     !print *, m
    end do
    !call cpu_time(finish1)
    !print '("Time = ",f6.3," seconds.")',finish1-start1
    creal1(:,p_c)     = ONE/avg!SUM(x33, DIM=2)!ONE
    creal2(:,p_c)     = TWO/avg!SUM(x44, DIM=2)!TWO
    creal3(:,p_c)     = THREE/avg!SUM(x99, DIM=2)!THREE
    creal4(:,p_c)     = FOUR/avg!SUM(x100, DIM=2)!FOUR
    !pollkg1(:,p_c)    = FIVE/avg!SUM(fx1y1, DIM=2)!FIVE
    !pollkg2(:,p_c)    = SIX/avg!SUM(fx2y2, DIM=2)!SIX
    creal8(:,p_c)     = SEVEN/avg!SUM(fx1x2, DIM=2)!SEVEN
    creal9(:,p_c)     = EIGHT/avg!SUM(fy1y2, DIM=2)!EIGHT
    !creal5(:,p_c)     = NINE/avg!SUM(fcross1, DIM=2)!NINE
    !creal6(:,p_c)     = TEN/avg!SUM(fcross2, DIM=2)!TEN
    
    !print *, creal1(:,p_c)!avg!ONE
    

    !ONE   =       (/ (0, I = 1,fftsize) /)!0
    !TWO   =       (/ (0, I = 1,fftsize) /)!0
    !THREE =       (/ (0, I = 1,fftsize) /)!0
    !FOUR  =       (/ (0, I = 1,fftsize) /)!0
    !!FIVE  =       (/ (0, I = 1,fftsize) /)!(0,0)
    !!SIX   =       (/ (0, I = 1,fftsize) /)!(0, 0)
    !SEVEN =       (/ (0, I = 1,fftsize) /)!(0,0)
    !EIGHT =       (/ (0, I = 1,fftsize) /)!(0,0)



    ONE   =       0.0
    TWO   =       0.0
    THREE =       0.0
    FOUR  =       0.0
    !FIVE  =       (0,0)
    !SIX   =       (0,0)
    SEVEN =       (0.0,0.0)
    EIGHT =       (0.0,0.0)
    !NINE  =       (0.0,0.0)
    !TEN   =       (0.0,0.0)
    !call cpu_time(finish1)
    cpu_t= cpu_t + finish1-start1
    !print '("Time = ",f6.3," seconds.")',finish1-start1
    !print *, 'Time for FFT1 is..', cpu_t
  end do
  !print *, p_c 
  !print *, 'Time for FFT1 is..', cpu_t
  !creal8     =creal8/sqrt(creal1*creal3)
  !creal9     =creal9/sqrt(creal2*creal4)

  !creal5     =creal5/sqrt(creal1*creal4)
  !creal6     =creal6/sqrt(creal2*creal3)

  print *, 'Now destroying all the plans..'
  call cpu_time(finish)
  call dfftw_destroy_plan_ ( plan_forward1 )
  print '("Time = ",f6.3," seconds.")',finish-start

 
  end subroutine external_loop


