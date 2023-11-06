  !
  PROGRAM raman
  !
  USE io_raman, ONLY : stdout
  USE raman_mod
  USE raman_read
  USE raman_output
  !USE elopt
  USE calc_raman_intensity
  USE knum
  !
  IMPLICIT NONE
  !
  ! For date and time
  integer :: ti(1:8), i
  character(len=10) :: t(1:3)
  !
  ! initial message
  WRITE(stdout, '(/5x,a)') repeat('=',67)
  WRITE(stdout,'(a)')"                                                                            "
  WRITE(stdout,'(a)')"  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  "
  WRITE(stdout,'(a)')"  RRRRRRRRR         A        MM            MM         A        NNN      NN  "
  WRITE(stdout,'(a)')"  RR      RR       AAA       MMM          MMM        AAA       NNNN     NN  "
  WRITE(stdout,'(a)')"  RR      RR      AA AA      MMMM        MMMM       AA AA      NN NN    NN  "
  WRITE(stdout,'(a)')"  RRRRRRRRR      AA   AA     MM MM      MM MM      AA   AA     NN  NN   NN  "
  WRITE(stdout,'(a)')"  RRRR          AAAAAAAAA    MM  MM    MM  MM     AAAAAAAAA    NN   NN  NN  "
  WRITE(stdout,'(a)')"  RR  RR       AA       AA   MM   MM  MM   MM    AA       AA   NN    NN NN  "
  WRITE(stdout,'(a)')"  RR    RR    AA         AA  MM    MMMM    MM   AA         AA  NN     NNNN  "
  WRITE(stdout,'(a)')"  RR      RR AA           AA MM     MM     MM  AA           AA NN      NNN  "
  WRITE(stdout,'(a)')"  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  "
  WRITE(stdout,'(a)')"                                                                            "
  WRITE(stdout,'(5x,"Welcome to the Raman program: QERaman")') 
  WRITE(stdout,'(5x,"Please cite our paper at https://doi.org/10.1016/j.cpc.2023.108967")')
  WRITE(stdout,'(/5x,a)') repeat('=',67)
  !
  CALL date_and_time(t(1),t(2),t(3),ti)
  WRITE(stdout,'(A,i0,A,i0,A,i0,A,i0,A,i0,A,i0,A)') "     == Program started at    ", &
       ti(1),"/",ti(2),"/",ti(3),"  ",ti(5),":",ti(6),":",ti(7)," ==" 
  !
  ! Initial setting for parameters
  CALL init_parameter
  !
  ! Reading the input file
  CALL raman_readin
  !CALL raman_readdata
  !
  ! Initial setting of polarization vectors
  CALL init_polvec
  !
  ! Reading of matrix element
  CALL read_dvec
  CALL read_elph
  !
  ! Standard output of reading data
  CALL raman_out
  !
  ! Interpolate the kpoint
  IF (n_kinterp .GT. 0) THEN
     DO i = 1,n_kinterp
        CALL k_interpolate
     END DO
  END IF
  !
  CALL k_weight
  !
  ! Calculation of Raman intensity
  CALL raman_intensity
  !
  ! Output of the data file
  CALL raman_out_data
  !
  CALL date_and_time(t(1),t(2),t(3),ti)
  WRITE(stdout, '(1x)')
  WRITE(stdout,'(A,i0,A,i0,A,i0,A,i0,A,i0,A,i0,A)') "     == Program finished at    ", &
       ti(1),"/",ti(2),"/",ti(3),"  ",ti(5),":",ti(6),":",ti(7)," =="
  !
  WRITE(stdout, '(/5x,a)') repeat('=',67)
  WRITE(stdout, '(5x,"JOB DONE")')
  WRITE(stdout, '(/5x,a)') repeat('=',67)

  END PROGRAM raman
