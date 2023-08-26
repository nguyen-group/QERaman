  !
  ! Started by Y. Tatsumi 2016/12/22
  !
  MODULE knum
  ! Module for the interpolate and increase the k point
  !
  USE io_raman, ONLY : stdout
  USE raman_mod, ONLY : outdir, sorb, circular_pol, &
       nks, nbnd, nqs, nmode, nrs, nel, npol, &
       k, q, wk, eq, rs, elaser, eigv, intensity_raman, &
       dvec, matele_opt, matele_elph, &
       pi, im, &
       ev2cm, ry2cm, ry2ev
  !
  IMPLICIT NONE
  CONTAINS
  !-------------------------------------------------------------------
  SUBROUTINE k_interpolate
  !-------------------------------------------------------------------
  !
  ! Interpolation of k point
  !
  INTEGER :: ik, ik2, ip, ibi, ibf, i, ixyz
  INTEGER :: nk_interpolate, nks2, nks3, nky
  REAL(KIND=KIND(1.0d0)) :: eps=1.0d-6, iprod
  REAL(KIND=KIND(1.0d0)), ALLOCATABLE :: k2(:,:), k3(:,:), eigv2(:,:), eigv3(:,:), wk2(:)
  COMPLEX(KIND=KIND(1.0d0)), ALLOCATABLE :: dvec2(:,:,:,:), matele_elph2(:,:,:,:,:), dvec3(:,:,:,:), matele_elph3(:,:,:,:,:)
  CHARACTER(len=256) :: filename_kpoint
  LOGICAL :: line1
  !
  WRITE(stdout,'(1x)')
  WRITE(stdout,'(5x,"k points interpolation")')
  WRITE(stdout,'(1x)')
  !
  DO  ik = 1, nks
     DO ibi = 1, nbnd
        DO ibf = 1, nbnd
           !DO ixyz = 1, 3
              IF (REAL(dvec(ik,ibi,ibf,2)) .LT. 0.0d0) THEN
                 dvec(ik,ibi,ibf,:) = -dvec(ik,ibi,ibf,:)
              END IF
              !IF (AIMAG(dvec(ik,ibi,ibf,2)) .LT. 0.0d0) THEN
              !   dvec(ik,ibi,ibf,:) = CONJG(dvec(ik,ibi,ibf,:))
              !END IF              
           !END DO
        END DO
     END DO
  END DO
  ALLOCATE(k2(2*nks,3),dvec2(2*nks,nbnd,nbnd,3),matele_elph2(nmode,nqs,2*nks,nbnd,nbnd),eigv2(2*nks,nbnd))
  k2(:,:) = 0.0d0
  dvec2(:,:,:,:) = 0.0d0
  matele_elph2(:,:,:,:,:) = 0.0d0
  eigv2(:,:) = 0.0d0
  nks2 = 0
  DO ik = 1, nks
     IF (ik .EQ. 1) THEN
        nks2 = nks2 + 1
        k2(1,:) = k(1,:)
        dvec2(1,:,:,:) = dvec(1,:,:,:)
        matele_elph2(:,:,1,:,:) = matele_elph(:,:,1,:,:)
        eigv2(1,:) = eigv(1,:)
     ELSE IF (k(ik,1) .GT. k(ik-1,1)+eps) THEN
        nks2 = nks2+1
        k2(nks2,:) = k(ik,:)
        dvec2(nks2,:,:,:) = dvec(ik,:,:,:)
        matele_elph2(:,:,nks2,:,:) = matele_elph(:,:,ik,:,:)
        eigv2(nks2,:) = eigv(ik,:)
     ELSE
        nks2 = nks2 + 1
        k2(nks2,:) = (k(ik,:)+k(ik-1,:))/2.0d0
        dvec2(nks2,:,:,:) = ((dvec(ik,:,:,:)+dvec(ik-1,:,:,:))/2.0d0) * ((ABS(dvec(ik,:,:,:))+&
           ABS(dvec(ik-1,:,:,:)))/2.0d0) / (ABS((dvec(ik,:,:,:)+dvec(ik-1,:,:,:))/2.0d0))
        matele_elph2(:,:,nks2,:,:) = (matele_elph(:,:,ik,:,:)+matele_elph(:,:,ik-1,:,:))/2.0d0
        eigv2(nks2,:) = (eigv(ik,:)+eigv(ik-1,:))/2.0d0
        nks2 = nks2 + 1
        k2(nks2,:) = k(ik,:)
        dvec2(nks2,:,:,:) = dvec(ik,:,:,:)
        matele_elph2(:,:,nks2,:,:) = matele_elph(:,:,ik,:,:)
        eigv2(nks2,:) = eigv(ik,:)
     END IF
  END DO
  ALLOCATE(k3(nks2,3),dvec3(nks2,nbnd,nbnd,3),matele_elph3(nmode,nqs,nks2,nbnd,nbnd),eigv3(nks2,nbnd))
  k3(:,:) = 0.0d0
  dvec3(:,:,:,:) = 0.0d0
  matele_elph3(:,:,:,:,:) = 0.0d0
  eigv3(:,:) = 0.0d0
  DO ik = 1, nks2
     k3(ik,:) = k2(ik,:)
     dvec3(ik,:,:,:) = dvec2(ik,:,:,:)
     matele_elph3(:,:,ik,:,:) = matele_elph2(:,:,ik,:,:)
     eigv3(ik,:) = eigv2(ik,:)
  END DO
  DEALLOCATE(k2,dvec2,matele_elph2,eigv2)

  ALLOCATE(k2(2*nks2,3),dvec2(2*nks2,nbnd,nbnd,3),matele_elph2(nmode,nqs,2*nks2,nbnd,nbnd),eigv2(2*nks2,nbnd))
  k2(:,:) = 0.0d0
  dvec2(:,:,:,:) = 0.0d0
  matele_elph2(:,:,:,:,:) = 0.0d0
  eigv2(:,:) = 0.0d0
  nks3 = 0
  line1 = .TRUE.
  nky = 0
  DO ik = 1, nks2
     IF (k3(ik,1) .GT. k3(ik-1,1)+eps .AND. line1 .EQV. .TRUE.) THEN
        IF (nky .EQ. 0) THEN
           nky = ik-1
        END IF
        line1 = .FALSE.
        !write(*,*) nky
     END IF
     IF (line1 .EQV. .TRUE.) THEN
        nks3 = nks3+1
        k2(nks3,:) = k3(ik,:)
        dvec2(nks3,:,:,:) = dvec3(ik,:,:,:)
        matele_elph2(:,:,nks3,:,:) = matele_elph3(:,:,ik,:,:)
        eigv2(nks3,:) = eigv3(ik,:)
     ELSE IF (line1 .EQV. .FALSE.) THEN
        DO ik2 = ik, ik+nky-1
           nks3 = nks3 + 1
           k2(nks3,:) = (k3(ik2,:)+k3(ik2-nky,:))/2.0d0
           dvec2(nks3,:,:,:) = ((dvec3(ik2,:,:,:)+dvec3(ik2-nky,:,:,:))/2.0d0) * ((ABS(dvec3(ik2,:,:,:))+&
              ABS(dvec3(ik2-nky,:,:,:)))/2.0d0) / (ABS((dvec3(ik2,:,:,:)+dvec3(ik2-nky,:,:,:))/2.0d0))
           matele_elph2(:,:,nks3,:,:) = (matele_elph3(:,:,ik2,:,:)+matele_elph3(:,:,ik2-nky,:,:))/2.0d0
           eigv2(nks3,:) = (eigv3(ik2,:)+eigv3(ik2-nky,:))/2.0d0
        END DO
        nks3 = nks3+1
        k2(nks3,:) = k3(ik,:)
        dvec2(nks3,:,:,:) = dvec3(ik,:,:,:)
        matele_elph2(:,:,nks3,:,:) = matele_elph3(:,:,ik,:,:)
        eigv2(nks3,:) = eigv3(ik,:)
        line1 = .TRUE.
     END IF
  END DO
  !stop
  DEALLOCATE(k3,dvec3,matele_elph3,eigv3)
  ALLOCATE(k3(nks3,3),dvec3(nks3,nbnd,nbnd,3),matele_elph3(nmode,nqs,nks3,nbnd,nbnd),eigv3(nks3,nbnd))
  k3(:,:) = 0.0d0
  dvec3(:,:,:,:) = 0.0d0
  matele_elph3(:,:,:,:,:) = 0.0d0
  eigv3(:,:) = 0.0d0
  DO ik = 1, nks3
     k3(ik,:) = k2(ik,:)
     dvec3(ik,:,:,:) = dvec2(ik,:,:,:)
     matele_elph3(:,:,ik,:,:) = matele_elph2(:,:,ik,:,:)
     eigv3(ik,:) = eigv2(ik,:)
  END DO
  DEALLOCATE(k2,dvec2,matele_elph2,eigv2)
  
  !ku = 2.0d0 * pi / alat
  !filename_kpoint = 'kpoint_test.dat'
  !OPEN(16,file=TRIM(ADJUSTL(filename_kpoint)))
  !WRITE(16,'("#K point list (unit:2pi/a)")')
  !WRITE(16,'("#ik,kx,ky,kz")')
  !DO ik = 1, nks3
  !   WRITE(16,'(1x,I6,1x,3(1x,f10.4))') ik, k3(ik,1)/ku, k3(ik,2)/ku, k3(ik,3)/ku
  !ENDDO
  !CLOSE(16)

  !STOP
  !
  DEALLOCATE(k,dvec,matele_elph,eigv)
  nks = nks3
  ALLOCATE(k(nks,3),dvec(nks,nbnd,nbnd,3),matele_elph(nmode,nqs,nks,nbnd,nbnd),eigv(nks,nbnd))
  k(:,:) = 0.0d0
  dvec(:,:,:,:) = 0.0d0
  matele_elph(:,:,:,:,:) = 0.0d0
  eigv(:,:) = 0.0d0
  DO ik = 1, nks
     k(ik,:) = k3(ik,:)
     dvec(ik,:,:,:) = dvec3(ik,:,:,:)
     matele_elph(:,:,ik,:,:) = matele_elph3(:,:,ik,:,:)
     eigv(ik,:) = eigv3(ik,:)
  END DO
  ALLOCATE(wk2(nks))
  wk2(:) = wk(1)
  DEALLOCATE(wk)
  ALLOCATE(wk(nks))
  wk(:) = wk2(:)
  !
  DEALLOCATE(k3,dvec3,matele_elph3,eigv3,wk2)
END SUBROUTINE k_interpolate
!
!
SUBROUTINE k_weight
!
  INTEGER :: ik, ik2, ik3
  INTEGER :: nky
  REAL(KIND=KIND(1.0d0)) :: eps=1.0d-6
!
  nky = 0
  DO ik = 1, nks
     IF (ik .GT. 1 ) THEN
        IF(k(ik,1) .GT. k(ik-1,1)) THEN
           ik2 = ik
           EXIT
        END IF
     END IF
     wk(ik) = wk(ik)/2.0d0
     nky = nky + 1
  END DO
  !write(*,*) nky
  DO ik = ik2, nks
     IF(k(ik,1) .GT. k(ik-1,1)) THEN
        wk(ik) = wk(ik)/2.0d0
     END IF
     IF(k(ik+1,1) .GT. k(ik,1)) THEN
        wk(ik) = wk(ik)/2.0d0
     END IF
     IF(ik+nky .GT. nks) THEN
        ik3 = ik
        EXIT
     END IF
  END DO
  !write(*,*) nks
  !write(*,*) ik3
  DO ik = ik3,nks
     wk(ik) = wk(ik)/2.0d0
  END DO
  wk(1) = wk(1)/2.0d0
  wk(nky) = wk(nky)/2.0d0
  wk(nks-nky+1) = wk(nks-nky+1)/2.0d0
  wk(nks) = wk(nks)/2.0d0
  !write(*,*) wk(:)
!  
END SUBROUTINE k_weight
!
!
END MODULE knum
