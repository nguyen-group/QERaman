  !
  !This module is used for numbering the output in the program.
  !Created by Tatsymi 2016/11/15
  !Modifed by N. T. Hung 2020/10/10
  !-------------------------------------------------------------------
  MODULE io_raman
  !-------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  PRIVATE
  SAVE
  !
  PUBLIC :: stdout, dvecin, elphin, &
       raman_spec, kpoint, fn_matele_opt, fn_matele_elph, fn_raman_k, fn_eigv
  !
  ! Output of standard output (6-)
  !
  INTEGER :: stdout = 6 ! Standard output
  !
  ! Reading data file (101-)
  !
  !
  INTEGER :: dvecin = 101 ! data file of matrix element of dipole vector
  !
  INTEGER :: elphin = 102 ! data file of electron-phonon interaction
  !
  !
  ! Output of data files (201-)
  !
  INTEGER :: raman_spec = 201 ! Output for Raman spectra
  !
  INTEGER :: kpoint = 202 ! Output for calculated k point list
  !
  INTEGER :: fn_matele_opt = 203 ! Output for the electron-photon matrix element
  !
  INTEGER :: fn_matele_elph = 204 ! Output for the electron-phonon matrix element
  !
  INTEGER :: fn_raman_k = 205 ! Output for Raman matrix element
  !
  INTEGER :: fn_eigv = 206 ! Output for eigen energy list
  !
  END MODULE io_raman
