  !
  ! Started by Y. Tatsumi 2016/11/29
  ! Last modified by Y. Tatsumi 2017/09/05
  ! Last modified by N. T. Hung 2022/5/26
  !
  MODULE calc_raman_intensity
  ! Module to calculate the Raman intensity
  !
  USE  io_raman 
  USE  raman_mod, ONLY : nks, nbnd, nspin, npol, nbnd_occ, nmode, nqs, nel, nrs, &
       gamma, gamma_raman, rs_start, rs_end, &
       k, wk, eigv, eq, raman_k, intensity_raman, rs, &
       dvec, polvec, rtensor, rtensor_k, &
       matele_elph, elaser, &
       pi, im, &
       ry2ev, ev2ry, ev2cm, ry2cm
  !
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = KIND(1.d0)
  !
  CONTAINS
  !-------------------------------------------------------------------
  SUBROUTINE raman_intensity
    !-------------------------------------------------------------------
    !
    !  This subroutine calculate the Raman intensity.  
    !
    !
    INTEGER :: ibi, ibn1, ibn2, ibn3, ibf,ik, iq, imode, ipi, ips, irs, iel, i, j
    REAL(DP) :: r
    COMPLEX(DP), ALLOCATABLE :: matele_raman(:,:,:,:,:), rtensor_interv(:,:,:,:,:), rtensor_interv_k(:,:,:,:,:,:)
    !, rtensorb(:,:,:,:),matele_ramanb(:,:,:,:)
    !
    !
    ! for testing
    !matele_elph(:,:,:,:,:) = 1.0d0
    !dvec(:,:,:,:) = 1.0d0
    !matele_opt(:,:,:,:) = 1.0d0
    !
    ALLOCATE(matele_raman(nel,npol,npol,nmode,nqs), raman_k(nel,npol,npol,nmode,nks,nqs))
    matele_raman(:,:,:,:,:) = 0.0d0
    raman_k(:,:,:,:,:,:) = 0.0d0
    ALLOCATE(intensity_raman(nel,npol,npol,nrs))
    intensity_raman(:,:,:,:) = 0.0d0
    ALLOCATE(rs(nrs))
    rs(:) = 0.0d0
    !ALLOCATE(rtensor(nel,nmode,nqs,3,3),rtensor_k(nel,nmode,nks,nqs,3,3),rtensorb(nel,nqs,3,3),matele_ramanb(nel,npol,npol,nqs))
    ALLOCATE(rtensor(nel,nmode,nqs,3,3),rtensor_k(nel,nmode,nks,nqs,3,3),rtensor_interv(nel,nmode,nqs,3,3),rtensor_interv_k(nel,nmode,nks,nqs,3,3))
    rtensor(:,:,:,:,:) = 0.0d0
    rtensor_interv(:,:,:,:,:) = 0.0d0
    rtensor_interv_k(:,:,:,:,:,:) = 0.0d0
    rtensor_k(:,:,:,:,:,:) = 0.0d0
    !rtensorb(:,:,:,:) = 0.0d0
    !matele_ramanb(:,:,:,:) = 0.0d0
    !
    !DO ipi = 1, npol
    !   DO ips = 1, npol
    DO iel = 1, nel
       DO imode = 1, nmode
          DO iq = 1, nqs
             DO ik = 1, nks
                DO ibi = 1, nbnd_occ !install and final state
                   !DO ibi = nbv,nbv
                   DO ibn1 = nbnd_occ+1, nbnd !medium state
                      !DO ibn1 = nbv+1,nbv+1
                      DO ibn2 = nbnd_occ+1, nbnd !visual medium state
                         !DO ibn2 = nbv+1,nbv+1   
                         CALL raman_tensor(iel,imode,iq,ik,ibi,ibn1,ibn2)
                         !DO ibn3 = nbv+1,nbnd
                         !DO ibf = 1,nbv
                         !call raman_tensor_interv(iel,imode,iq,ik,ibi,ibn1,ibn2,ibn3,ibf)
                         !END DO
                         !END DO
                      END DO ! Loop for ibn2
                   END DO ! Loop for ibn1
                END DO ! Loop for ibi
                !DO i = 1, 3
                !   DO j = 1, 3
                rtensor_k(iel,imode,ik,iq,:,:) = rtensor_k(iel,imode,ik,iq,:,:) * wk(ik)
                rtensor(iel,imode,iq,:,:) = rtensor(iel,imode,iq,:,:) + rtensor_k(iel,imode,ik,iq,:,:)
                !rtensorb(iel,iq,:,:) = rtensorb(iel,iq,:,:) + rtensor_k(iel,imode,ik,iq,:,:)
                !   END DO
                !END DO
             END DO ! Loop for ik
          END DO ! Loop for iq
       END DO ! Loop for imode
    END DO ! Loop for iel
    !  END DO ! Loop for ips
    !END DO !Loop for ipi
    !
    DO ipi = 1, npol
       DO ips = 1, npol
          !      DO iel = 1, nel
          !         DO imode = 1, nmode
          !            DO iq = 1, nqs
          DO i = 1, 3
             DO j = 1, 3
                !                     DO ik = 1, nks
                raman_k(:,ipi,ips,:,:,:) = raman_k(:,ipi,ips,:,:,:) + CONJG(polvec(ips,i))*rtensor_k(:,:,:,:,i,j)*polvec(ipi,j)
                !END DO
                matele_raman(:,ipi,ips,:,:) = matele_raman(:,ipi,ips,:,:) + CONJG(polvec(ips,i))*rtensor(:,:,:,i,j)*polvec(ipi,j)
                !matele_ramanb(:,ipi,ips,:) =  matele_ramanb(:,ipi,ips,:) + CONJG(polvec(ips,i))*rtensorb(:,:,i,j)*polvec(ipi,j)
             END DO
          END DO
          !           END DO
          !        END DO
          !     END DO
       END DO
    END DO
    !raman_k(:,:,:,:,:,:) = raman_k(:,:,:,:,:,:) / DBLE(nks)
    !matele_raman(:,:,:,:,:) = matele_raman(:,:,:,:,:) / DBLE(nks)
    !matele_ramanb(:,:,:,:) = matele_ramanb(:,:,:,:) / DBLE(nks)
    DEALLOCATE(rtensor_k)
    !
    DO iel = 1, nel
       DO imode = 1, nmode
          DO i = 1, 3
             DO j = 1, 3
                IF (ABS(rtensor(iel,imode,1,i,j)) .GT. 1.0d0) THEN
                   rtensor(iel,imode,1,:,:) = rtensor(iel,imode,1,:,:) / ABS(rtensor(iel,imode,1,i,j))
                END IF
             END DO
          END DO
       END DO
    END DO
    !
    ! Standard output of Raman tensor
    WRITE(stdout, '(1x)')
    WRITE(stdout, '(1x)') 
    WRITE(stdout, '(5x,"Raman tensor")')
    !
    DO iel = 1, nel
       WRITE(stdout, '(5x,A6,3x,I2)') "iel = ", iel
       DO imode = 1, nmode
          WRITE(stdout, '(7x,A6,3x,I2)') "imode = ", imode
          WRITE(stdout, '(9x,6(2x,e16.8e3))') rtensor(iel,imode,1,1,1), rtensor(iel,imode,1,1,2), rtensor(iel,imode,1,1,3)
          WRITE(stdout, '(9x,6(2x,e16.8e3))') rtensor(iel,imode,1,2,1), rtensor(iel,imode,1,2,2), rtensor(iel,imode,1,2,3)
          WRITE(stdout, '(9x,6(2x,e16.8e3))') rtensor(iel,imode,1,3,1), rtensor(iel,imode,1,3,2), rtensor(iel,imode,1,3,3)
       END DO
    END DO
    WRITE(stdout, '(1x)')
    WRITE(stdout, '(1x)')
    DEALLOCATE(rtensor)
    !
    DO ipi = 1, npol
       DO ips = 1, npol       
          DO iel = 1, nel
             DO iq = 1, nqs
                DO irs = 1, nrs
                   rs(irs) = rs_start+ (DBLE(irs-1)/DBLE(nrs-1))*(rs_end-rs_start)
                   !
                   r = 0.d0
                   DO imode = 1, nmode
                      !r = r + matele_raman(iel,ipi,ips,imode,iq) * (1.0d0/(sqrt(2.0d0*pi)*gamma_raman)) * exp(-((eq(imode,iq)-rs(irs))**2.0d0) / (2.0d0 * (gamma_raman)**2.0d0))
                      !rs(irs) = 2.0d-4 * DBLE(irs)
                      !
                      r = r + ((DBLE(CONJG(matele_raman(iel,ipi,ips,imode,iq))*matele_raman(iel,ipi,ips,imode,iq))) &
                      * (1.0d0/pi) * (gamma_raman/((eq(imode,iq)-rs(irs))**2.0d0+gamma_raman**2.0d0))) 
                      !
                      ! by Lorentzian function (see page 2 in PRB 97, 115407 (2018))
                      ! intensity_raman(iel,ipi,ips,irs) = intensity_raman(iel,ipi,ips,irs) + r &
                      !  * (1.0d0/(sqrt(2.0d0*pi)*gamma_raman)) * exp(-((eq(imode,iq)-rs(irs))**2.0d0) / (2.0d0 * (gamma_raman)**2.0d0)) ) 
                      ! by Gaussian
                      !
                      ! take sum_k|Raman|^2 for testing
                      !DO ik = 1, nks
                      !  r = r + (DBLE(CONJG(raman_k(iel,ipi,ips,imode,ik,iq))*raman_k(iel,ipi,ips,imode,ik,iq)) &
                      !     * (1.0d0/pi) * (gamma_raman/((eq(imode,iq)-rs(irs))**2.0d0+gamma_raman**2.0d0)))
                      !ENDDO
                   END DO ! Loop for imode
                   intensity_raman(iel,ipi,ips,irs) = r
                END DO ! Loop for irs
             END DO ! Loop for iq
          END DO ! Loop for iel
       END DO ! Loop for ips
    END DO ! Loop for ipi
    DEALLOCATE(eq,matele_raman)
    !
  END SUBROUTINE raman_intensity
  !
  !
  !-------------------------------------------------------------------
  SUBROUTINE raman_tensor(iel,imode,iq,ik,ibi,ibn1,ibn2)
  !-------------------------------------------------------------------
  !
  ! This subroutine construct Raman tensor
  !
  INTEGER :: iel, imode, iq, ik, ibi, ibn1, ibn2
  !
  INTEGER :: i, j
  COMPLEX(DP) :: rt
  !
  !
  !
  DO i = 1, 3
     DO j = 1, 3
        ! Product of three matrix element
        ! [D(k,f2,i)*M_q(v,k,f1,f2)*D(k,i,f1)]
        ! Raman intensity should proportion to 1/E_L^4 since M_opt proportion to D/E_L
        rt = (dvec(ik,ibn2,ibi,i) / (elaser(iel)-eq(imode,iq))) &
             * matele_elph(imode,iq,ik,ibn1,ibn2) * (dvec(ik,ibi,ibn1,j) / elaser(iel))
        !rt = dvec(ik,ibn2,ibi,i) * matele_elph(imode,iq,ik,ibn1,ibn2) * dvec(ik,ibi,ibn1,j)
        ! Non-resonant
        !  rt = rt / ((eigv(ik,ibn1)-eigv(ik,ibi)) * (eigv(ik,ibn2)-eigv(ik,ibi)))
        ! Resonant condition
        rt = rt / ((elaser(iel)-(eigv(ik,ibn1)-eigv(ik,ibi))-im*gamma) &
             * (elaser(iel)-(eigv(ik,ibn2)-eigv(ik,ibi))-eq(imode,iq)-im*gamma))
        !make Raman tensor
        rtensor_k(iel,imode,ik,iq,i,j) = rtensor_k(iel,imode,ik,iq,i,j) + rt
     END DO
  END DO
  !
  !           
  END SUBROUTINE raman_tensor
  !
  !
  !
  END MODULE calc_raman_intensity
    
  
