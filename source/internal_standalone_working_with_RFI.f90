subroutine external_loop_rfi (comf_X, comf_Y, comf1_X, comf1_Y, avg, le, fftsize, RFI_X1, RFI_Y1, RFI_X2, RFI_Y2, len_1, len_2, &
 & creal1, creal2, creal3, creal4, creal5, creal6, creal8, creal9, pollkg1, pollkg2)

! Computes Correlations.
!
! Author:
!   Pavan Uttarkar  
  
  integer(kind=4), intent(in):: len_1
  integer(kind=4), intent(in):: len_2
  integer(kind=4), intent(in):: fftsize
  integer(kind=4), intent(in):: avg
  integer(kind=4), intent(in):: le
  
  integer(kind=1), intent(in), dimension(0:len_1)::comf_X
  integer(kind=1), intent(in), dimension(0:len_1)::comf_Y
  integer(kind=1), intent(in), dimension(0:len_2)::comf1_X
  integer(kind=1), intent(in), dimension(0:len_2)::comf1_Y

  !RFI Matrix!
  integer(kind=1), intent(in), dimension(0:fftsize)::RFI_X1
  integer(kind=1), intent(in), dimension(0:fftsize)::RFI_Y1
  integer(kind=1), intent(in), dimension(0:fftsize)::RFI_X2
  integer(kind=1), intent(in), dimension(0:fftsize)::RFI_Y2


  !Auto Correlation!
  real             , intent(out), dimension(0:fftsize, 0:le)::creal1
  real             , intent(out), dimension(0:fftsize, 0:le)::creal2
  real             , intent(out), dimension(0:fftsize, 0:le)::creal3
  real             , intent(out), dimension(0:fftsize, 0:le)::creal4

  !Cross Pol Correlation!
  double complex, intent(out), dimension(0:fftsize, 0:le)::creal5
  double complex, intent(out), dimension(0:fftsize, 0:le)::creal6

  !Cross Pol leakage!
  double complex, intent(out), dimension(0:fftsize, 0:le)::pollkg1
  double complex, intent(out), dimension(0:fftsize, 0:le)::pollkg2

  !Cross Correlation!
  double complex, intent(out), dimension(0:fftsize, 0:le)::creal8
  double complex, intent(out), dimension(0:fftsize, 0:le)::creal9
 
  integer(kind=4) :: p_c
  integer(kind=4)  ::avg1
  !Intermediate products..!
  real          , dimension(0:fftsize, 0:avg    ):: x33
  real          , dimension(0:fftsize, 0:avg    ):: x44
  real          , dimension(0:fftsize, 0:avg    ):: x99
  real          , dimension(0:fftsize, 0:avg    ):: x100
  double complex, dimension(0:fftsize, 0:avg    ):: fx1y1
  double complex, dimension(0:fftsize, 0:avg    ):: fx2y2
  double complex, dimension(0:fftsize, 0:avg    ):: fx1x2
  double complex, dimension(0:fftsize, 0:avg    ):: fy1y2
  double complex, dimension(0:fftsize, 0:avg    ):: fcross1
  double complex, dimension(0:fftsize, 0:avg    ):: fcross2

  integer(kind=8)  ::fir
  integer(kind=8)  ::sed
  integer(kind=8)  ::sed1
  integer(kind=8)  ::m
  integer(kind=8)  ::le1

  integer(kind=1), dimension(0:fftsize*avg*2)::chunk_X
  integer(kind=1), dimension(0:fftsize*avg*2):: chunk_Y
  integer(kind=1), dimension(0:fftsize*avg*2):: chunk1_X
  integer(kind=1), dimension(0:fftsize*avg*2):: chunk1_Y

  real   (kind=8), dimension(0:fftsize*2)  ::X1 
  real   (kind=8), dimension(0:fftsize*2)  ::Y1
  real   (kind=8), dimension(0:fftsize*2)  ::X2
  real   (kind=8), dimension(0:fftsize*2)  ::Y2


  double complex, dimension(0:fftsize*2)  ::fX1
  double complex, dimension(0:fftsize*2)  ::fY1
  double complex, dimension(0:fftsize*2)  ::fX2
  double complex, dimension(0:fftsize*2)  ::fY2

  double complex, dimension(0:fftsize)  ::fX1X1
  double complex, dimension(0:fftsize)  ::fY1Y1
  double complex, dimension(0:fftsize)  ::fX2X2
  double complex, dimension(0:fftsize)  ::fY2Y2


  double complex, dimension(0:fftsize)  ::fX1_1
  double complex, dimension(0:fftsize)  ::fY1_1
  double complex, dimension(0:fftsize)  ::fX2_1
  double complex, dimension(0:fftsize)  ::fY2_1

  real(kind=8)  , dimension(0:fftsize)  ::ONE  
  real(kind=8)  , dimension(0:fftsize)  ::TWO  
  real(kind=8)  , dimension(0:fftsize)  ::THREE
  real(kind=8)  , dimension(0:fftsize)  ::FOUR 
  double complex, dimension(0:fftsize)  ::FIVE 
  double complex, dimension(0:fftsize)  ::SIX  
  double complex, dimension(0:fftsize)  ::SEVEN
  double complex, dimension(0:fftsize)  ::EIGHT
  double complex, dimension(0:fftsize)  ::NINE 
  double complex, dimension(0:fftsize)  ::TEN  

  integer(kind=8)::plan_forward1, plan_forward2, plan_forward3, plan_forward4
  !integer(kind=1)::iret
  !call dfftw_init_threads(iret)
  !print *, iret
  !call dfftw_plan_with_nthreads(8)
  sed1  =  avg*fftsize*2
  avg1    = avg-1


  le1 = le-1
    
  print *, 'Planning FFT'
  call dfftw_plan_dft_r2c_1d ( plan_forward1, 512, X1, fX1, FFTW_ESTIMATE )!call dfftw_plan_dft_r2c_1d_ ( plan_forward1, 512, X1, fX1, FFTW_ESTIMATE )
  call dfftw_plan_dft_r2c_1d ( plan_forward2, 512, Y1, fY1, FFTW_ESTIMATE )!call dfftw_plan_dft_r2c_1d_ ( plan_forward2, 512, Y1, fY1, FFTW_ESTIMATE )
  call dfftw_plan_dft_r2c_1d ( plan_forward3, 512, X2, fX2, FFTW_ESTIMATE )!call dfftw_plan_dft_r2c_1d_ ( plan_forward3, 512, X2, fX2, FFTW_ESTIMATE )
  call dfftw_plan_dft_r2c_1d ( plan_forward4, 512, Y2, fY2, FFTW_ESTIMATE )!call dfftw_plan_dft_r2c_1d_ ( plan_forward4, 512, Y2, fY2, FFTW_ESTIMATE )
  ONE   =       0
  TWO   =       0
  THREE =       0
  FOUR  =       0
  FIVE  =       (0,0)
  SIX   =       (0, 0)
  SEVEN =       (0,0)
  EIGHT =       (0,0)
  NINE  =       (0,0)
  TEN   =       (0,0)

  do p_c =0, le1
    !print *, p_c
    fir             = avg*p_c*fftsize*2
    sed             = avg*(p_c+1)*fftsize*2

    chunk_X         = comf_X(fir:sed)
    chunk1_X        = comf1_X(fir:sed)
    chunk_Y         = comf_Y(fir:sed)
    chunk1_Y        = comf1_Y(fir:sed)


    do m=0, avg1
      !write(, '(A)', advance='no') m!print *, !'In secondary loop..'
      X1  =  chunk_X(m*fftsize*2:(m+1)*fftsize*2)
      Y1  =  chunk_Y(m*fftsize*2:(m+1)*fftsize*2)
      X2  =  chunk1_X(m*fftsize*2:(m+1)*fftsize*2)
      Y2  =  chunk1_Y(m*fftsize*2:(m+1)*fftsize*2)
  
      call dfftw_execute( plan_forward1 )!call dfftw_execute_ ( plan_forward1 )
      fX1_1 = fX1(0:255)
      fX1_1 = fX1_1*RFI_X1
      fX1X1    = fX1_1*CONJG(fX1_1)
      x33(:,m) = fX1X1
      

      call dfftw_execute( plan_forward2 ) !call dfftw_execute_ ( plan_forward2 )
      fY1_1 = fY1(0:255)
      fY1_1 = fY1_1*RFI_Y1
      fY1Y1    = fY1_1*CONJG(fY1_1)
      x44(:,m) = fY1Y1

      call dfftw_execute( plan_forward3 ) !call dfftw_execute_ ( plan_forward3 )
      fX2_1 = fX2(0:255)
      fX2_1 = fX2_1*RFI_X2
      fX2X2    = fX2_1*CONJG(fX2_1)
      x99(:,m) = fX2X2

      call dfftw_execute( plan_forward4 ) !call dfftw_execute_ ( plan_forward4 )
      fY2_1 = fY2(0:255)
      fX2_1 = fX2_1*RFI_Y2
      fY2Y2     = fY2_1*CONJG(fY2_1)
      x100(:,m) = fY2Y2
      
      fx1x2(:,m)  =  fX1_1*conjg(fX2_1)
      fy1y2(:,m)  =  fY1_1*conjg(fY2_1)

      fcross1(:,m) =  conjg(fX1_1)*fY2_1
      fcross2(:,m) =  conjg(fX2_1)*fY1_1

      fx1y1(:,m)   = conjg(fX1_1)*fY1_1
      fx2y2(:,m)   = conjg(fX2_1)*fY2_1

      ONE          = ONE   + x33(:,m)
      TWO          = TWO   + x44(:,m)
      THREE        = THREE + x99(:,m)
      FOUR         = FOUR  + x100(:,m)
      FIVE         = FIVE  + fx1y1(:,m)
      SIX          = SIX   + fx2y2(:,m)
      SEVEN        = SEVEN + fx1x2(:,m)
      EIGHT        = EIGHT + fy1y2(:,m)
      NINE         = NINE  + fcross1(:,m)
      TEN          = TEN   + fcross2(:,m)

    end do
    
    creal1(:,p_c)     = ONE
    creal2(:,p_c)     = TWO
    creal3(:,p_c)     = THREE
    creal4(:,p_c)     = FOUR
    pollkg1(:,p_c)    = FIVE
    pollkg2(:,p_c)    = SIX
    creal8(:,p_c)     = SEVEN
    creal9(:,p_c)     = EIGHT
    creal5(:,p_c)     = NINE
    creal6(:,p_c)     = TEN
    
    ONE   =       0
    TWO   =       0
    THREE =       0
    FOUR  =       0
    FIVE  =       (0,0)
    SIX   =       (0,0)
    SEVEN =       (0,0)
    EIGHT =       (0,0)
    NINE  =       (0,0)
    TEN   =       (0,0)

  end do

  !print *, p_c 
  creal8     =creal8/sqrt(creal1*creal3)
  creal9     =creal9/sqrt(creal2*creal4)

  creal5     =creal5/sqrt(creal1*creal4)
  creal6     =creal6/sqrt(creal2*creal3)

  print *, 'Now destroying all the plans..'
  !call dfftw_destroy_plan_ ( plan_forward1 )
  !call dfftw_destroy_plan_ ( plan_forward2 )
  !call dfftw_destroy_plan_ ( plan_forward3 )
  !call dfftw_destroy_plan_ ( plan_forward4 )
 
  end subroutine external_loop_rfi


