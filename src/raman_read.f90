  !
  ! Started by Y. Tatsumi 2016/11/15
  ! Modifed by N. T. Hung 2022/10/10
  MODULE raman_read
  ! Module to read input and data files
  !
  USE  io_raman, ONLY : stdout, dvecin, elphin
  USE  raman_mod, ONLY : prefix, outdir, fil_dvec, fil_elph, &
       sorb, circular_pol, nonpol, plot_matele_opt, plot_matele_elph, plot_raman_k, &
       nks, nksq, nbnd, nspin, nmode, nqs, nel, nrs, nbnd_occ, n_kinterp, &
       gamma, gamma_raman, rs_start, rs_end, &
       k, q, wk, wq, eigv, eq, &
       matele_elph, dvec, elaser, &
       pi, im, ev2ry, ry2cm
  !
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = KIND(1.d0)
  !
  CONTAINS
  !-------------------------------------------------------------------
  SUBROUTINE raman_readin
  !-------------------------------------------------------------------
  !
  !  This subroutine reads the variable from the input file.  
  !
  !
  IMPLICIT NONE
  !
  INTEGER :: iel, ios
  REAL(DP) :: elaser1,elaser2,elaser3,elaser4,elaser5,elaser6,elaser7
  !
  !
  namelist / inputraman / &
       outdir, prefix, fil_dvec, fil_elph, &
       sorb, circular_pol, nonpol, plot_matele_opt, plot_matele_elph, plot_raman_k, &
       gamma, gamma_raman, rs_start, rs_end, nrs, &
       elaser1, elaser2, elaser3, elaser4, elaser5, elaser6, elaser7, n_kinterp
  !
  ! Set default values for variables in namelist
  !
  outdir = '.'
  prefix = ''
  fil_dvec = ''
  fil_elph = ''
  sorb = .FALSE.
  circular_pol = .FALSE.
  nonpol = .FALSE.
  plot_matele_opt = .FALSE.
  plot_matele_elph = .FALSE.
  plot_raman_k = .FALSE.
  rs_start = 0.0
  rs_end = 800.0
  nrs = 1000
  elaser1 = 2.33d0 !or 532 nm
  elaser2 = -1.0d0
  elaser3 = -1.0d0
  elaser4 = -1.0d0
  elaser5 = -1.0d0
  elaser6 = -1.0d0
  elaser7 = -1.0d0
  n_kinterp = 0
  !
  ! Read the namelist inputraman
  READ (5, inputraman, iostat = ios)
  !
  ! modify the unit to atomic unit
  elaser1 = elaser1 * ev2ry
  elaser2 = elaser2 * ev2ry
  elaser3 = elaser3 * ev2ry
  elaser4 = elaser4 * ev2ry
  elaser5 = elaser5 * ev2ry
  elaser6 = elaser6 * ev2ry
  elaser7 = elaser7 * ev2ry
  gamma = gamma * ev2ry
  gamma_raman = gamma_raman * ev2ry
  rs_start = rs_start / ry2cm
  rs_end = rs_end / ry2cm
  !
  nel = 0
  IF (elaser1 .GT. 0.0d0) THEN
     nel = nel + 1
  END IF
  IF (elaser2 .GT. 0.0d0) THEN
     nel = nel + 1
  END IF
  IF (elaser3 .GT. 0.0d0) THEN
     nel = nel + 1
  END IF
  IF (elaser4 .GT. 0.0d0) THEN
     nel = nel + 1
  END IF
  IF (elaser5 .GT. 0.0d0) THEN
     nel = nel + 1
  END IF
  IF (elaser6 .GT. 0.0d0) THEN
     nel = nel + 1
  END IF
  IF (elaser7 .GT. 0.0d0) THEN
     nel = nel + 1
  END IF
  ALLOCATE(elaser(nel))
  iel = 1
  IF (elaser1 .GT. 0.0d0) THEN
     elaser(iel) = elaser1
     iel = iel+1
  END IF
  IF (elaser2 .GT. 0.0d0) THEN
     elaser(iel) = elaser2
     iel = iel+1
  END IF
  IF (elaser3 .GT. 0.0d0) THEN
     elaser(iel) = elaser3
     iel = iel+1
  END IF
  IF (elaser4 .GT. 0.0d0) THEN
     elaser(iel) = elaser4
     iel = iel+1
  END IF
  IF (elaser5 .GT. 0.0d0) THEN
     elaser(iel) = elaser5
     iel = iel+1
  END IF
  IF (elaser6 .GT. 0.0d0) THEN
     elaser(iel) = elaser6
     iel = iel+1
  END IF
  IF (elaser7 .GT. 0.0d0) THEN
     elaser(iel) = elaser7
  END IF
  !
  !
  END SUBROUTINE raman_readin
  !
  !
  !-------------------------------------------------------------------
  SUBROUTINE read_dvec
  !-------------------------------------------------------------------
  !
  !  This subroutine reads the data file of matrix element of 
  !  dipole vector from modified bands.x in the QERaman code. 
  !
  ! Local prameters
  INTEGER :: ik, ii, jj, ribi, ribf, rik
  CHARACTER(LEN=256) :: dummy
  REAL(DP) :: rdx_r, rdx_i, rdy_r, rdy_i, rdz_r, rdz_i 
  !
  WRITE(stdout,'(1x)')
  WRITE(stdout,'(5x,"Reading dipole vector matrix element data file")')
  WRITE(stdout,'(1x)')
  !
  OPEN(dvecin, file=TRIM(ADJUSTL(fil_dvec)), FORM = 'formatted', &
        STATUS = 'old', ACTION = 'read', ACCESS = 'sequential')
  !
  READ(dvecin,*) dummy
  READ(dvecin,*) dummy
  READ(dvecin,*) nbnd, nbnd_occ, nks
  READ(dvecin,*) dummy
  READ(dvecin,*) dummy
  READ(dvecin,*) dummy
  !
  ALLOCATE(k(nks,3), wk(nks), eigv(nks,nbnd))
  ALLOCATE(dvec(nks,nbnd,nbnd,3))
  !
  DO ik = 1, nks
     !
     READ(dvecin, *) k(ik,:), wk(ik)
     !
     DO ii = 1,nbnd !final state
        DO jj = 1,nbnd ! install state
           READ(dvecin, *) ribi, ribf, eigv(ik,ii), eigv(ik,jj), rdx_r, rdx_i, rdy_r, rdy_i, rdz_r, rdz_i 
           dvec(ik,ii,jj,1) = rdx_r + im * rdx_i  !D_x^(k,f,i)
           dvec(ik,ii,jj,2) = rdy_r + im * rdy_i  !D_y^(k,f,i)
           dvec(ik,ii,jj,3) = rdz_r + im * rdz_i  !D_z^(k,f,i)
        ENDDO
     ENDDO
  ENDDO
  eigv(:,:) = eigv(:,:) * ev2ry
  !
  CLOSE(dvecin, STATUS = 'keep')
  !
  WRITE(stdout,'(1x)')
  WRITE(stdout,'(5x,"Finished reading dipole vector matrix element")')
  WRITE(stdout,'(1x)')
  !
  !
  END SUBROUTINE read_dvec
  !
  !
  !-------------------------------------------------------------------
  SUBROUTINE read_elph
  !-------------------------------------------------------------------
  !
  !  This subroutine reads the data file of electron-phonon interaction  
  !  from modified ph.x in the QERaman code. 
  !
  IMPLICIT NONE
  !
  ! Local prameters
  INTEGER :: riq, rimode, ribi, ribf, rik, rikq
  INTEGER :: iq, ik, ii, jj, nu
  REAL(DP) :: rkx, rky, rkz, rwk
  REAL(DP) :: reigvi, reigvj, rmatele, rmatele_r, rmatele_i
  CHARACTER(LEN=256) :: dummy
  !
  WRITE(stdout,'(1x)')
  WRITE(stdout,'(5x,"Reading electron-phonon matrix element data file")')
  WRITE(stdout,'(1x)')
  !
  OPEN(elphin, file = TRIM(ADJUSTL(fil_elph)), FORM = 'formatted', &
        STATUS = 'old', ACTION = 'read', ACCESS = 'sequential')
  !
  READ(elphin,*) dummy
  READ(elphin,*) dummy
  READ(elphin,*) nbnd, nmode, nqs, nks, nksq
  !write(*,*)nbnd, nmode, nqs, nks, nksq
  READ(elphin,*) dummy
  READ(elphin,*) dummy
  READ(elphin,*) dummy
  READ(elphin,*) dummy
  !
  ALLOCATE(q(nqs,3), wq(nqs), eq(nmode,nqs))
  !ALLOCATE(k(nks,3), wk(nks), eigv(nksq,nbnd))
  ALLOCATE(matele_elph(nmode,nqs,nksq,nbnd,nbnd))
  !
  DO iq = 1, nqs
     READ (elphin, *) q(iq,:), wq(iq), riq
     !write(*,*) q(iq,:), wq(iq), riq
     DO ik = 1, nksq
        !READ (elphin, *) k(ik,:), wq(ik), rik, rikq
        READ (elphin, *) rkx, rky, rkz, rwk, rik, rikq
        DO ii = 1, nbnd
           DO jj = 1, nbnd
              DO nu = 1, nmode
                 !READ (elphin, *) ribi, ribf, rimode, eigv(ik,ii), eigv(ik,jj), eq(nu,iq), rmatele, rmatele_r, rmatele_i
                 READ (elphin, *) ribi, ribf, rimode, reigvi, reigvj, eq(nu,iq), rmatele, rmatele_r, rmatele_i
                 matele_elph(nu,iq,ik,ii,jj) = rmatele_r + im * rmatele_i !M_q(v,k,i,f)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  ! Convert omega meV to Ry
  !eigv(:,:) = eigv(:,:) * ev2ry
  eq(:,:) = eq(:,:) * 1.0d-3 * ev2ry
  matele_elph(:,:,:,:,:) = matele_elph(:,:,:,:,:) * 1.0d-3 * ev2ry
  !
  CLOSE(elphin, STATUS = 'keep')
  !
  WRITE(stdout,'(1x)')
  WRITE(stdout,'(5x,"Finished reading electron-phonon matrix element")')
  WRITE(stdout,'(1x)')
  !
  END SUBROUTINE read_elph
  !
  !
  END MODULE raman_read


  
