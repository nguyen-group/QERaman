  !
  ! Started by Y. Tatsumi 2016/11/28
  ! Modified by N. T. Hung 2022/10/10
  MODULE raman_output
  ! Module for the output
  !
  USE io_raman, ONLY : stdout, raman_spec, kpoint, fn_eigv, fn_matele_opt, fn_matele_elph, fn_raman_k
  USE raman_mod, ONLY : outdir, sorb, circular_pol, nonpol, plot_matele_opt, plot_matele_elph, plot_raman_k, &
       nks, nbnd, nqs, nmode, nrs, nel, npol, nbnd_occ, &
       k, q, wk, wq, eigv, eq, rs, elaser, raman_k, intensity_raman, &
       dvec,matele_opt, matele_elph, polvec, &
       pi, &
       ev2cm, ry2cm, ry2ev
  !
  IMPLICIT NONE
  CONTAINS
  !-------------------------------------------------------------------
  SUBROUTINE raman_out
  !-------------------------------------------------------------------
  !
  ! Standard output for Raman program
  !
  INTEGER :: ik, iq, imode
  ! 
  WRITE(stdout, '(1x)')
  WRITE(stdout, '(5x,"Number of k points",I6)') nks
  WRITE(stdout, '(5x,"Number of q points",I6)') nqs
  WRITE(stdout, '(5x,"Number of bands",I6)') nbnd
  WRITE(stdout, '(5x,"Number of valence bands",I6)') nbnd_occ
  !WRITE(stdout, '(5x,"Number of G vectors",I6)') ngs
  WRITE(stdout, '(5x,"Number of phonon mode",I6)') nmode
  WRITE(stdout, '(1x)')
  !
  !WRITE(stdout, '(5x,"Lattice constant alat = ",f8.4,3x,"(au)")') alat
  !
  !WRITE(stdout, '(1x)')
  !
  !WRITE(stdout,'(5x,"Primitive lattice vectors (au)")')
  !WRITE(stdout,'(6x,A6,3(1x,f8.4),A1)') 'a1 = (', avec(1,1), avec(1,2), avec(1,3), ')'
  !WRITE(stdout,'(6x,A6,3(1x,f8.4),A1)') 'a2 = (', avec(2,1), avec(2,2), avec(2,3), ')'
  !WRITE(stdout,'(6x,A6,3(1x,f8.4),A1)') 'a3 = (', avec(3,1), avec(3,2), avec(3,3), ')'
  !WRITE(stdout,'(5x,"Reciprocal lattice vectors (au)")')
  !WRITE(stdout,'(6x,A6,3(1x,f8.4),A1)') 'b1 = (', bvec(1,1), bvec(1,2), bvec(1,3), ')'
  !WRITE(stdout,'(6x,A6,3(1x,f8.4),A1)') 'b2 = (', bvec(2,1), bvec(2,2), bvec(2,3), ')'
  !WRITE(stdout,'(6x,A6,3(1x,f8.4),A1)') 'b3 = (', bvec(3,1), bvec(3,2), bvec(3,3), ')'
  !
  !WRITE(stdout, '(1x)')
  !
  WRITE(stdout,'(5x,"k points (cart. coord. in units 2pi/alat)")')
  DO ik = 1, nks
     WRITE(stdout,'(6x,I5,3x,A1,3(1x,f8.4),A12,1x,f8.4)') ik,'(', k(ik,1), k(ik,2), k(ik,3), '),  wk = ', wk(ik)
  ENDDO
  !
  WRITE(stdout, '(1x)')
  !
  WRITE(stdout,'(5x,"q points and phonon eigen values (cm-1)")')
  DO iq = 1, nqs
     WRITE(stdout,'(6x,I5,3x,A1,3(1x,f8.4),A1)') iq,'(', q(iq,1), q(iq,2), q(iq,3), ')'
     WRITE(stdout,'(5x,1000f10.4)') (eq(imode,iq)*ry2cm,imode=1,nmode)
  ENDDO
  !
  !
  END SUBROUTINE raman_out
  !
  !
  !-------------------------------------------------------------------
  SUBROUTINE raman_out_data
  !-------------------------------------------------------------------
  !
  ! Output for the data files
  !  !
  INTEGER :: iel, imode, ik, irs, ibi, ibf, ip
  REAL(KIND=KIND(1.0d0)) :: eps=1.0d-6
  !
  CHARACTER(len=256) :: make_outdir, make_outdir_matele_opt, make_outdir_matele_elph, make_outdir_raman_k
  CHARACTER(len=256) :: outdir_matele_opt, outdir_matele_elph, outdir_raman_k
  CHARACTER(len=256) :: filename_raman_spec
  CHARACTER(len=256) :: filename_kpoint, filename_eigv
  CHARACTER(len=256) :: filename_matele_opt
  CHARACTER(len=256) :: filename_matele_elph
  CHARACTER(len=256) :: cel, cmopt1, cmopt2, cmelph1, cmelph2
  CHARACTER(len=256) :: filename_raman_k
  CHARACTER(len=256) :: cramank1, cramank2, eramank, mramank
  !
  WRITE(make_outdir,'("mkdir -p ",A228)') outdir
  !
  CALL SYSTEM(make_outdir)
  !
  ! Data output for Raman spectra
  DO iel = 1, nel
     WRITE(cel,'(I1)') iel
     filename_raman_spec = TRIM(ADJUSTL(outdir))//'/raman_spectra'//TRIM(ADJUSTL(cel))//'.dat'
     OPEN(raman_spec,file=TRIM(ADJUSTL(filename_raman_spec)))
     WRITE(raman_spec,'("#Raman spectra")')
     WRITE(raman_spec,'(A13,1x,f6.2,1x,A2)') '#Laser energy', elaser(iel)*ry2ev, 'eV'
     IF (circular_pol .EQV. .FALSE. .AND.  nonpol .EQV. .FALSE.) THEN
        WRITE(raman_spec,'("#Raman shift (cm-1), Raman Spectra (p=XX,YY,ZZ,XY,YX)")')
        DO irs = 1,nrs
           WRITE(raman_spec,'(1x,f10.4,5(1x,e16.8e3))') rs(irs)*ry2cm,intensity_raman(iel,1,1,irs), &
                intensity_raman(iel,2,2,irs),intensity_raman(iel,3,3,irs),intensity_raman(iel,1,2,irs), &
                intensity_raman(iel,2,1,irs)
        ENDDO
     ELSE IF (circular_pol .EQV. .FALSE. .AND.  nonpol .EQV. .TRUE.) THEN
        WRITE(raman_spec,'("#Raman shift (cm-1), Raman Spectra (p=XX,YY,ZZ,XY,YX,nonpol)")')
        DO irs = 1,nrs
           WRITE(raman_spec,'(1x,f10.4,6(1x,e16.8e3))') rs(irs)*ry2cm,intensity_raman(iel,1,1,irs), &
                intensity_raman(iel,2,2,irs),intensity_raman(iel,3,3,irs),intensity_raman(iel,1,2,irs), &
                intensity_raman(iel,2,1,irs),intensity_raman(iel,6,6,irs)
        ENDDO
     ELSE IF (circular_pol .EQV. .TRUE. .AND. nonpol .EQV. .FALSE.) THEN
        WRITE(raman_spec,'("#Raman shift (cm-1), Raman Spectra (XX,YY,ZZ,s+s+,s-s-,XY,YX,s+s-,s-s+)")')
        DO irs = 1,nrs
           WRITE(raman_spec,'(1x,f10.4,9(1x,e16.8e3))') rs(irs)*ry2cm,intensity_raman(iel,1,1,irs), &
                intensity_raman(iel,2,2,irs),intensity_raman(iel,3,3,irs),intensity_raman(iel,4,4,irs), &
                intensity_raman(iel,5,5,irs),intensity_raman(iel,1,2,irs),intensity_raman(iel,2,1,irs), &
                intensity_raman(iel,4,5,irs),intensity_raman(iel,5,4,irs)
        ENDDO
     ELSE IF (circular_pol .EQV. .TRUE. .AND. nonpol .EQV. .TRUE.) THEN
                WRITE(raman_spec,'("#Raman shift (cm-1), Raman Spectra (XX,YY,ZZ,s+s+,s-s-,XY,YX,s+s-,s-s+,nonpol)")')
        DO irs = 1,nrs
           WRITE(raman_spec,'(1x,f10.4,10(1x,e16.8e3))') rs(irs)*ry2cm,intensity_raman(iel,1,1,irs), &
                intensity_raman(iel,2,2,irs),intensity_raman(iel,3,3,irs),intensity_raman(iel,4,4,irs), &
                intensity_raman(iel,5,5,irs),intensity_raman(iel,1,2,irs),intensity_raman(iel,2,1,irs), &
                intensity_raman(iel,4,5,irs),intensity_raman(iel,5,4,irs),intensity_raman(iel,6,6,irs)
        ENDDO
     ENDIF
     CLOSE(raman_spec)
  ENDDO
  !
  ! Data output for kpoint list
  !ku = 2.0d0 * pi / alat
  filename_kpoint = TRIM(ADJUSTL(outdir))//'/kpoint.dat'
  OPEN(kpoint,file=TRIM(ADJUSTL(filename_kpoint)))
  !WRITE(kpoint,'("#K point list (unit:2pi/a)")')
  WRITE(kpoint,'("#K point list (cart. coord. in units 2pi/alat)")')
  WRITE(kpoint,'("#ik, kx, ky, kz, weight of k")')
  DO ik = 1, nks
     WRITE(kpoint,'(1x,I6,1x,4(1x,f10.4))') ik, k(ik,1), k(ik,2), k(ik,3), wk(ik)
  ENDDO
  CLOSE(kpoint)
  !
  ! Data output for eigen energy list
  filename_eigv = TRIM(ADJUSTL(outdir))//'/eigv.dat'
  OPEN(fn_eigv,file=TRIM(ADJUSTL(filename_eigv)))
  WRITE(fn_eigv,'("#Eigen energy list (unit:eV)")')
  WRITE(fn_eigv,'("#ik, kx, ky, kz, eigen energy(1:nbnd)")')
  DO ik = 1, nks
     !IF (ik > 1 .AND. (k(ik,1) > k(ik-1,1)+eps .OR. k(ik,1) < k(ik-1,1)-eps)) THEN
     !   WRITE(fn_eigv,*) ' '
     !END IF
     WRITE(fn_eigv,'(1x,I6,1x,3(1x,f10.4),1000(1x,e16.8e3))') ik, k(ik,1), k(ik,2), k (ik,3), eigv(ik,:)*ry2ev
  ENDDO
  CLOSE(kpoint)
  !
  ! Data output for electron-photon matrix element
  IF (plot_matele_opt .EQV. .TRUE.) THEN
     ALLOCATE(matele_opt(nks,nbnd,nbnd,npol))
     matele_opt(:,:,:,:) = 0.0d0
     outdir_matele_opt = TRIM(ADJUSTL(outdir))//'/matele_opt/'
     WRITE(make_outdir_matele_opt,'("mkdir -p ",A228)') outdir_matele_opt
     CALL SYSTEM(make_outdir_matele_opt)
     DO ibi = 1, nbnd
        WRITE(cmopt1,'(I3)') ibi
        DO ibf = 1, nbnd
           WRITE(cmopt2,'(I3)') ibf
           filename_matele_opt = TRIM(ADJUSTL(outdir_matele_opt))//'matele_opt_'//&
              &TRIM(ADJUSTL(cmopt1))//'_'//TRIM(ADJUSTL(cmopt2))//'.dat'
           OPEN(fn_matele_opt, file=filename_matele_opt)
           DO ik = 1, nks
              DO ip = 1, npol
                 matele_opt(ik,ibi,ibf,ip) = polvec(ip,1)*CONJG(dvec(ik,ibi,ibf,1)) &
                    + polvec(ip,2)*CONJG(dvec(ik,ibi,ibf,2)) + polvec(ip,3)*CONJG(dvec(ik,ibi,ibf,3))
              END DO
              !IF (ik > 1 .AND. (k(ik,1) > k(ik-1,1)+eps .OR. k(ik,1) < k(ik-1,1)-eps)) THEN
              !   WRITE(fn_matele_opt,*) ' '
              !END IF
              !write(301,'(3(1x,f10.4),5(1x,e16.8e3))') k(ik,1), k(ik,2), k (ik,3), abs(matele_opt(ik,4,5,1)), abs(matele_opt(ik,5,4,2)), abs(matele_opt(ik,5,4,3)), abs(matele_elph(5,1,ik,5,5)), abs(matele_elph(6,1,ik,5,5))
                   IF (circular_pol .EQV. .FALSE. .AND.  nonpol .EQV. .FALSE.) THEN
                      WRITE(fn_matele_opt,'(3(1x,f10.4),9(1x,e16.8e3))') k(ik,1), k(ik,2), k(ik,3), &
                         REAL(matele_opt(ik,ibi,ibf,1)), AIMAG(matele_opt(ik,ibi,ibf,1)), &
                         ABS(matele_opt(ik,ibi,ibf,1)), REAL(matele_opt(ik,ibi,ibf,2)), &
                         AIMAG(matele_opt(ik,ibi,ibf,2)), ABS(matele_opt(ik,ibi,ibf,2)), &
                         REAL(matele_opt(ik,ibi,ibf,3)), AIMAG(matele_opt(ik,ibi,ibf,3)), ABS(matele_opt(ik,ibi,ibf,3))
                   ELSE IF (circular_pol .EQV. .TRUE. .AND.  nonpol .EQV. .FALSE.) THEN
                      WRITE(fn_matele_opt,'(3(1x,f10.4),15(1x,e16.8e3))') k(ik,1), k(ik,2), k(ik,3), &
                         REAL(matele_opt(ik,ibi,ibf,1)), AIMAG(matele_opt(ik,ibi,ibf,1)), &
                         ABS(matele_opt(ik,ibi,ibf,1)), REAL(matele_opt(ik,ibi,ibf,2)), &
                         AIMAG(matele_opt(ik,ibi,ibf,2)), ABS(matele_opt(ik,ibi,ibf,2)), &
                         REAL(matele_opt(ik,ibi,ibf,3)), AIMAG(matele_opt(ik,ibi,ibf,3)), &
                         ABS(matele_opt(ik,ibi,ibf,3)), REAL(matele_opt(ik,ibi,ibf,4)), &
                         AIMAG(matele_opt(ik,ibi,ibf,4)), ABS(matele_opt(ik,ibi,ibf,4)), &
                         REAL(matele_opt(ik,ibi,ibf,5)), AIMAG(matele_opt(ik,ibi,ibf,5)), ABS(matele_opt(ik,ibi,ibf,5))
                   ELSE
                      !WRITE(fn_matele_opt,'(3(1x,f10.4),18(1x,e16.8e3))') k(ik,1), k(ik,2), k(ik,3), REAL(matele_opt(ik,ibi,ibf,1)), AIMAG(matele_opt(ik,ibi,ibf,1)), ABS(matele_opt(ik,ibi,ibf,1)), REAL(matele_opt(ik,ibi,ibf,2)), AIMAG(matele_opt(ik,ibi,ibf,2)), ABS(matele_opt(ik,ibi,ibf,2)), REAL(matele_opt(ik,ibi,ibf,3)), AIMAG(matele_opt(ik,ibi,ibf,3)), ABS(matele_opt(ik,ibi,ibf,3)), REAL(matele_opt(ik,ibi,ibf,4)), AIMAG(matele_opt(ik,ibi,ibf,4)), ABS(matele_opt(ik,ibi,ibf,4)), REAL(matele_opt(ik,ibi,ibf,5)), AIMAG(matele_opt(ik,ibi,ibf,5)), ABS(matele_opt(ik,ibi,ibf,5)), REAL(matele_opt(ik,ibi,ibf,6)), AIMAG(matele_opt(ik,ibi,ibf,6)), ABS(matele_opt(ik,ibi,ibf,6))
                      WRITE(fn_matele_opt,'(3(1x,f10.4),22(1x,e16.8e3))') k(ik,1), k(ik,2), k(ik,3), &
                         REAL(matele_opt(ik,ibi,ibf,1)), AIMAG(matele_opt(ik,ibi,ibf,1)), &
                         ABS(matele_opt(ik,ibi,ibf,1)), REAL(matele_opt(ik,ibi,ibf,2)), &
                         AIMAG(matele_opt(ik,ibi,ibf,2)), ABS(matele_opt(ik,ibi,ibf,2)), &
                         REAL(matele_opt(ik,ibi,ibf,3)), AIMAG(matele_opt(ik,ibi,ibf,3)), &
                         ABS(matele_opt(ik,ibi,ibf,3)), REAL(matele_opt(ik,ibi,ibf,4)), &
                         AIMAG(matele_opt(ik,ibi,ibf,4)), ABS(matele_opt(ik,ibi,ibf,4)), &
                         REAL(matele_opt(ik,ibi,ibf,5)), AIMAG(matele_opt(ik,ibi,ibf,5)), &
                         ABS(matele_opt(ik,ibi,ibf,5)), REAL(matele_opt(ik,ibi,ibf,6)), &
                         AIMAG(matele_opt(ik,ibi,ibf,6)), ABS(matele_opt(ik,ibi,ibf,6)), &
                         REAL(matele_opt(ik,ibi,ibf,4)*CONJG(matele_opt(ik,ibi,ibf,5))), &
                         AIMAG(matele_opt(ik,ibi,ibf,4)*CONJG(matele_opt(ik,ibi,ibf,5))), &
                         REAL(matele_opt(ik,ibi,ibf,4)*CONJG(matele_opt(ik,ibi,ibf,4))), &
                         AIMAG(matele_opt(ik,ibi,ibf,4)*CONJG(matele_opt(ik,ibi,ibf,4)))
                   END IF
           END DO
        END DO
     END DO
  END IF
  !
  ! Data output for electron-phonon matrix element
  IF (plot_matele_elph .EQV. .TRUE.) THEN
     outdir_matele_elph = TRIM(ADJUSTL(outdir))//'/matele_elph/'
     WRITE(make_outdir_matele_elph,'("mkdir -p ",A228)') outdir_matele_elph
     CALL SYSTEM(make_outdir_matele_elph)
     DO ibi = 1, nbnd
        WRITE(cmelph1,'(I3)') ibi
        DO ibf = 1, nbnd
           WRITE(cmelph2,'(I3)') ibf
           filename_matele_elph = TRIM(ADJUSTL(outdir_matele_elph))//'matele_elph_'//&
              &TRIM(ADJUSTL(cmelph1))//'_'//TRIM(ADJUSTL(cmelph2))//'.dat'
           OPEN(fn_matele_elph, file=filename_matele_elph)
           DO ik = 1, nks
              !IF (ik > 1 .AND. (k(ik,1) > k(ik-1,1)+eps .OR. k(ik,1) < k(ik-1,1)-eps)) THEN
              !   WRITE(fn_matele_elph,*) ' '
              !END IF
              WRITE(fn_matele_elph,'(3(1x,f10.4),1000(1x,e16.8e3))') k(ik,1), k(ik,2), k (ik,3), &
                 abs(matele_elph(:,1,ik,ibi,ibf)), real(matele_elph(:,1,ik,ibi,ibf)), aimag(matele_elph(:,1,ik,ibi,ibf))
           END DO
        END DO
     END DO
  END IF
  !
  !
  ! Data output of Raman matrix element for each k point
  IF (plot_raman_k .EQV. .TRUE.) THEN
     outdir_raman_k = TRIM(ADJUSTL(outdir))//'/raman_k/'
     WRITE(make_outdir_raman_k,'("mkdir -p ",A228)') outdir_raman_k
     CALL SYSTEM(make_outdir_raman_k)
     DO iel = 1, nel
        WRITE(eramank,'(I3)') iel
        DO imode = 1, nmode
           WRITE(mramank,'(I3)') imode
           !DO ibi = 1, nbnd
              !WRITE(cramank1,'(I3)') ibi
              !DO ibf = 1, nbnd
                 !WRITE(cramank2,'(I3)') ibf
                 filename_raman_k = TRIM(ADJUSTL(outdir_raman_k))//'raman_k_mode'//&
                    &TRIM(ADJUSTL(mramank))//'_elaser'//TRIM(ADJUSTL(eramank))//'.dat'
                 OPEN(fn_raman_k, file=filename_raman_k)
                 DO ik = 1, nks
                    !IF (ik > 1 .AND. (k(ik,1) > k(ik-1,1)+eps .OR. k(ik,1) < k(ik-1,1)-eps)) THEN
                    !   WRITE(fn_raman_k,*) ' '
                    !END IF
                    IF (circular_pol .EQV. .FALSE. .AND.  nonpol .EQV. .FALSE.) THEN
                       WRITE(fn_raman_k,'(3(1x,f10.4),10(1x,e16.8e3))')  k(ik,1), k(ik,2), k (ik,3), &
                            REAL(raman_k(iel,1,1,imode,ik,1)), AIMAG(raman_k(iel,1,1,imode,ik,1)), &
                            REAL(raman_k(iel,2,2,imode,ik,1)), AIMAG(raman_k(iel,2,2,imode,ik,1)), &
                            REAL(raman_k(iel,3,3,imode,ik,1)), AIMAG(raman_k(iel,3,3,imode,ik,1)), &
                            REAL(raman_k(iel,1,2,imode,ik,1)), AIMAG(raman_k(iel,1,2,imode,ik,1)), &
                            REAL(raman_k(iel,2,1,imode,ik,1)), AIMAG(raman_k(iel,2,1,imode,ik,1))
                    ELSE IF (circular_pol .EQV. .FALSE. .AND.  nonpol .EQV. .TRUE.) THEN
                       WRITE(fn_raman_k,'(3(1x,f10.4),12(1x,e16.8e3))')  k(ik,1), k(ik,2), k (ik,3), &
                            REAL(raman_k(iel,1,1,imode,ik,1)), AIMAG(raman_k(iel,1,1,imode,ik,1)), &
                            REAL(raman_k(iel,2,2,imode,ik,1)), AIMAG(raman_k(iel,2,2,imode,ik,1)), &
                            REAL(raman_k(iel,3,3,imode,ik,1)), AIMAG(raman_k(iel,3,3,imode,ik,1)), &
                            REAL(raman_k(iel,1,2,imode,ik,1)), AIMAG(raman_k(iel,1,2,imode,ik,1)), &
                            REAL(raman_k(iel,2,1,imode,ik,1)), AIMAG(raman_k(iel,2,1,imode,ik,1)), &
                            REAL(raman_k(iel,6,6,imode,ik,1)), AIMAG(raman_k(iel,6,6,imode,ik,1))
                    ELSE IF (circular_pol .EQV. .TRUE. .AND. nonpol .EQV. .FALSE.) THEN
                       WRITE(fn_raman_k,'(3(1x,f10.4),18(1x,e16.8e3))')  k(ik,1), k(ik,2), k (ik,3), &
                            REAL(raman_k(iel,1,1,imode,ik,1)), AIMAG(raman_k(iel,1,1,imode,ik,1)), &
                            REAL(raman_k(iel,2,2,imode,ik,1)), AIMAG(raman_k(iel,2,2,imode,ik,1)), &
                            REAL(raman_k(iel,3,3,imode,ik,1)), AIMAG(raman_k(iel,3,3,imode,ik,1)), &
                            REAL(raman_k(iel,4,4,imode,ik,1)), AIMAG(raman_k(iel,4,4,imode,ik,1)), &
                            REAL(raman_k(iel,5,5,imode,ik,1)), AIMAG(raman_k(iel,5,5,imode,ik,1)), &
                            REAL(raman_k(iel,1,2,imode,ik,1)), AIMAG(raman_k(iel,1,2,imode,ik,1)), &
                            REAL(raman_k(iel,2,1,imode,ik,1)), AIMAG(raman_k(iel,2,1,imode,ik,1)), &
                            REAL(raman_k(iel,4,5,imode,ik,1)), AIMAG(raman_k(iel,4,5,imode,ik,1)), &
                            REAL(raman_k(iel,5,4,imode,ik,1)), AIMAG(raman_k(iel,5,4,imode,ik,1))
                    ELSE IF (circular_pol .EQV. .TRUE. .AND. nonpol .EQV. .TRUE.) THEN
                       WRITE(fn_raman_k,'(3(1x,f10.4),20(1x,e16.8e3))')  k(ik,1), k(ik,2), k (ik,3), &
                            REAL(raman_k(iel,1,1,imode,ik,1)), AIMAG(raman_k(iel,1,1,imode,ik,1)), &
                            REAL(raman_k(iel,2,2,imode,ik,1)), AIMAG(raman_k(iel,2,2,imode,ik,1)), &
                            REAL(raman_k(iel,3,3,imode,ik,1)), AIMAG(raman_k(iel,3,3,imode,ik,1)), &
                            REAL(raman_k(iel,4,4,imode,ik,1)), AIMAG(raman_k(iel,4,4,imode,ik,1)), &
                            REAL(raman_k(iel,5,5,imode,ik,1)), AIMAG(raman_k(iel,5,5,imode,ik,1)), &
                            REAL(raman_k(iel,1,2,imode,ik,1)), AIMAG(raman_k(iel,1,2,imode,ik,1)), &
                            REAL(raman_k(iel,2,1,imode,ik,1)), AIMAG(raman_k(iel,2,1,imode,ik,1)), &
                            REAL(raman_k(iel,4,5,imode,ik,1)), AIMAG(raman_k(iel,4,5,imode,ik,1)), &
                            REAL(raman_k(iel,5,4,imode,ik,1)), AIMAG(raman_k(iel,5,4,imode,ik,1)), &
                            REAL(raman_k(iel,6,6,imode,ik,1)), AIMAG(raman_k(iel,6,6,imode,ik,1))
                    END IF
                 END DO
              !END DO
           !END DO
        END DO
     END DO
  END IF
  !
  !
  END SUBROUTINE raman_out_data
  !
  !
  !
  END MODULE raman_output
  
