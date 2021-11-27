subroutine external_loop (comf_X, comf_Y, avg, le, RFI1, RFI2, fftsize,  len_1,  &
 & creal1, creal2, creal3, creal4)
  
  integer(kind=8), intent(in):: len_1
  integer(kind=8), intent(in):: fftsize
  integer(kind=8), intent(in):: avg
  integer(kind=8), intent(in):: le
  
  integer(kind=1), intent(in), dimension(len_1)::comf_X
  integer(kind=1), intent(in), dimension(len_1)::comf_Y

  integer(kind=1), intent(in), dimension(fftsize)::RFI1
  integer(kind=1), intent(in), dimension(fftsize)::RFI2

  !Auto Correlation!
  real(kind=8)             , intent(out), dimension(fftsize, le)::creal1
  real(kind=8)             , intent(out), dimension(fftsize, le)::creal2
  real(kind=8)             , intent(out), dimension(fftsize, le)::creal3
  real(kind=8)             , intent(out), dimension(fftsize, le)::creal4


  integer(kind=8) :: p_c
  !Intermediate products..!
  real(kind=8)          , dimension(fftsize, avg    ):: x33
  real(kind=8)          , dimension(fftsize, avg    ):: x44


  integer(kind=8)  ::fir
  integer(kind=8)  ::sed
  integer(kind=8)  ::fir1
  integer(kind=8)  ::sed1
  integer(kind=8)  ::m

  integer(kind=1), dimension(fftsize*avg*2)::chunk_X
  integer(kind=1), dimension(fftsize*avg*2):: chunk_Y

  real   (kind=8), dimension(fftsize*2)  ::X1 


  double complex, dimension(fftsize*2)  ::fX1
  double complex, dimension(fftsize)  ::fX1X1
  double complex, dimension(fftsize)  ::fY1Y1


  double complex, dimension(fftsize)  ::fX1_1
  double complex, dimension(fftsize)  ::fY1_1



  integer(kind=8)::plan_forward1!, plan_forward2
  sed1  =  avg*fftsize*2
    
  print *, 'Planning FFT'
  call dfftw_plan_dft_r2c_1d ( plan_forward1, 512, X1, fX1, FFTW_ESTIMATE )!call dfftw_plan_dft_r2c_1d_ ( plan_forward1, 512, X1, fX1, FFTW_ESTIMATE )

  do p_c =1, le
    !print *, p_c
    fir             = avg*(p_c-1)*fftsize*2+1
    sed             = avg*(p_c)*fftsize*2
    !print*, fir, ":", sed
    chunk_X         = comf_X(fir:sed)
    chunk_Y         = comf_Y(fir:sed)


    do m=1, avg
      fir1=  (m-1)*fftsize*2+1
      sed1=  m*fftsize*2
      X1  =  chunk_X(fir1:sed1)



      call dfftw_execute( plan_forward1 )
      fX1_1       = fX1(1:256)
      fX1X1       = fX1_1*CONJG(fX1_1)
      x33(:,m)    = REAL( fX1X1 )



      X1  =  chunk_Y(fir1:sed1)
      call dfftw_execute( plan_forward1 )
      fY1_1 = fX1(1:256)
      fY1Y1 = fY1_1*CONJG(fY1_1)
      x44(:,m) = REAL( fY1Y1 )


    end do
    
    creal1(:,p_c)     = SUM(x33, DIM=2)
    creal1(:,p_c)     = creal1(:,p_c)*RFI1!ONE(1:256)
    creal2(:,p_c)     = SUM(x44, DIM=2)
    creal1(:,p_c)     = creal1(:,p_c)*RFI2!TWO(1:256)
    creal3(:,p_c)     = SUM(x33, DIM=2)
    creal4(:,p_c)     = SUM(x44, DIM=2)

  end do

  print *, 'Now destroying all the plans..'
  call dfftw_destroy_plan_ ( plan_forward1 )
  
end subroutine! external_loop
