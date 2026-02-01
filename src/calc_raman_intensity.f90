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
       dvec, polvec, rtensor, rtensor_k, rtensor_lpl, rtensor_llp, rtensor_pll, &
       matele_elph, elaser, temp, efermi, &
       kb, pi, im, &
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
    ALLOCATE(rtensor(nel,nmode,nqs,3,3), rtensor_k(nel,nmode,nks,nqs,3,3),&
       &rtensor_interv(nel,nmode,nqs,3,3), rtensor_interv_k(nel,nmode,nks,nqs,3,3))
    ALLOCATE(rtensor_lpl(nel,nmode,nks,nqs,3,3))
    ALLOCATE(rtensor_llp(nel,nmode,nks,nqs,3,3))
    ALLOCATE(rtensor_pll(nel,nmode,nks,nqs,3,3))
    rtensor(:,:,:,:,:) = 0.0d0
    rtensor_interv(:,:,:,:,:) = 0.0d0
    rtensor_interv_k(:,:,:,:,:,:) = 0.0d0
    rtensor_k(:,:,:,:,:,:) = 0.0d0
    rtensor_lpl(:,:,:,:,:,:) = 0.0d0
    rtensor_llp(:,:,:,:,:,:) = 0.0d0
    rtensor_pll(:,:,:,:,:,:) = 0.0d0
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
                rtensor_lpl(iel,imode,ik,iq,:,:) = rtensor_lpl(iel,imode,ik,iq,:,:) * wk(ik)
                rtensor_llp(iel,imode,ik,iq,:,:) = rtensor_llp(iel,imode,ik,iq,:,:) * wk(ik)
                rtensor_pll(iel,imode,ik,iq,:,:) = rtensor_pll(iel,imode,ik,iq,:,:) * wk(ik)
                !
                rtensor_k(iel,imode,ik,iq,:,:) = rtensor_lpl(iel,imode,ik,iq,:,:) &
                                                + rtensor_llp(iel,imode,ik,iq,:,:) &
                                                + rtensor_pll(iel,imode,ik,iq,:,:)
                !
                rtensor(iel,imode,iq,:,:) = rtensor(iel,imode,iq,:,:) + rtensor_k(iel,imode,ik,iq,:,:)
                !
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
    DEALLOCATE(rtensor_lpl)
    DEALLOCATE(rtensor_llp)
    DEALLOCATE(rtensor_pll)
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
  REAL(DP) :: nq, fi, fn1, fn2, chi, xi, fa
  !
  INTEGER :: i, j
  COMPLEX(DP) :: rt_lpl, rte_lpl, rth_lpl
  COMPLEX(DP) :: rt_llp, rte_llp, rth_llp
  COMPLEX(DP) :: rt_pll, rte_pll, rth_pll
  !
  ! Bose Einstein distribution for phonon and Fermi distribution for electrons
  ! avoid numerical error at T = 0 K
  IF (temp <= 1.0D-8) THEN
     nq = 0.0D0
     fi  = MERGE(1.0D0, 0.0D0, eigv(ik,ibi) < efermi)  !final state
     fn1 = MERGE(1.0D0, 0.0D0, eigv(ik,ibn1) < efermi) !intermediate state
     fn2 = MERGE(1.0D0, 0.0D0, eigv(ik,ibn2) < efermi) !scatterted state
  ELSE
     nq  = 1.0D0/(EXP(eq(imode,iq)/(kb*temp))-1.0D0)
     fi  = 1.0D0/(EXP((eigv(ik,ibi)-efermi)/(kb*temp))+1.0D0)
     fn1 = 1.0D0/(EXP((eigv(ik,ibn1)-efermi)/(kb*temp))+1.0D0)
     fn2 = 1.0D0/(EXP((eigv(ik,ibn2)-efermi)/(kb*temp))+1.0D0)
  END IF
  ! avoid acounting error at eq = 0
  IF (eq(imode,iq) < 1.0D-10) THEN
     nq = 0.0D0
  ELSE
     nq = 1.0D0/(EXP(eq(imode,iq)/(kb*temp))-1.0D0)
  END IF
  !
  ! Factor for electron-phonon scattering process [PRB 90, 035443 (2014)]
  IF (.NOT.((ibi == ibn1) .OR. (ibi == ibn2) .OR. (ibn1 == ibn2))) THEN
     chi = SQRT(nq + 1.D0) * fi * (1.D0 - fn1) * (1.D0 - fn2)
  ELSE
     chi = SQRT(nq + 1.D0) * fi * (1.D0 - fn2)
  END IF
  ! Factor for hole-phonon scattering process [PRB 90, 035443 (2014)]
  IF (.NOT.((ibi == ibn1) .OR. (ibi == ibn2) .OR. (ibn1 == ibn2))) THEN
     xi = SQRT(nq + 1.D0) * fi * fn1 * (1.D0 - fn2)
  ELSE
     xi = 0.d0
  END IF
  !
  ! Factor for Raman intensity calculation, in which Raman intensity is inversely proportional to the fourth power of the laser energy
  ! Ref. R. Loudon (1964) The Raman effect in crystals, Advances in Physics, 13:52, 423-482, DOI: 10.1080/00018736400101051
  fa = (elaser(iel)-eq(imode,iq))*elaser(iel)
  !
  DO i = 1, 3
     DO j = 1, 3
        !
        ! Product of Raman matrix element for photon absorption–phonon emission–photon creation (rpr)
        ! Add contribution from electron process [see PRB 90, 035443 (2014)]
        ! i*gamma means each energy transition must decay forward in time
        rte_lpl = dvec(ik,ibi,ibn2,i) * matele_elph(imode,iq,ik,ibn2,ibn1) * dvec(ik,ibn1,ibi,j) &
                / ((eigv(ik,ibn2)-eigv(ik,ibi)+im*gamma+eq(imode,iq)-elaser(iel))    &
                         * (eigv(ik,ibn1)-eigv(ik,ibi)+im*gamma-elaser(iel))*fa)
        ! Add contribution from hole process
        rth_lpl = dvec(ik,ibn1,ibn2,i) * matele_elph(imode,iq,ik,ibi,ibn1) * dvec(ik,ibn2,ibi,j) &
                / ((eigv(ik,ibn2)-eigv(ik,ibn1)+im*gamma+eq(imode,iq)-elaser(iel))    &
                         * (eigv(ik,ibn2)-eigv(ik,ibi)+im*gamma-elaser(iel))*fa)
        ! Total Raman tensor element for rpr process
        rt_lpl = chi*rte_lpl + xi*rth_lpl
        ! make Raman tensor for rpr process
        rtensor_lpl(iel,imode,ik,iq,i,j) = rtensor_lpl(iel,imode,ik,iq,i,j) + rt_lpl
        !
        ! Product of Raman matrix element for photon absorption–photon emission–phonon creation (rrp)
        ! Add contribution from electron process
        rte_llp = matele_elph(imode,iq,ik,ibi,ibn2) * dvec(ik,ibn2,ibn1,i) * dvec(ik,ibn1,ibi,j) &
                / ((eigv(ik,ibn1)-eigv(ik,ibi)+im*gamma-elaser(iel))          &
                         * (eigv(ik,ibn2)-eigv(ik,ibi)+im*gamma-eq(imode,iq))*fa)
        ! Add contribution from hole process
        rth_llp = matele_elph(imode,iq,ik,ibn1,ibn2) * dvec(ik,ibi,ibn1,i) * dvec(ik,ibn2,ibi,j) &
                / ((eigv(ik,ibn2)-eigv(ik,ibi)+im*gamma-elaser(iel))         &
                         * (eigv(ik,ibn2)-eigv(ik,ibn1)+im*gamma-eq(imode,iq))*fa)
        ! Total Raman tensor element for rrp process
        rt_llp = chi*rte_llp + xi*rth_llp
        ! make Raman tensor for rrp process
        rtensor_llp(iel,imode,ik,iq,i,j) = rtensor_llp(iel,imode,ik,iq,i,j) + rt_llp
        !
        ! Product of Raman matrix element for phonon creation–photon absorption–photon emission (prr)
        ! Add contribution from electron process
        rte_pll = dvec(ik,ibi,ibn2,i) * dvec(ik,ibn2,ibn1,j) * matele_elph(imode,iq,ik,ibn1,ibi) &
                / ((eigv(ik,ibn2)-eigv(ik,ibi)+im*gamma+eq(imode,iq)-elaser(iel))          &
                         * (eigv(ik,ibn1)-eigv(ik,ibi)+im*gamma+eq(imode,iq))*fa)
        ! Add contribution from hole process
        rth_pll = dvec(ik,ibn1,ibn2,i) * dvec(ik,ibi,ibn1,j) * matele_elph(imode,iq,ik,ibn2,ibi) &
                / ((eigv(ik,ibn2)-eigv(ik,ibn1)+im*gamma+eq(imode,iq)-elaser(iel))         &
                         * (eigv(ik,ibn2)-eigv(ik,ibi)+im*gamma+eq(imode,iq))*fa)
        ! Total Raman tensor element for prr process
        rt_pll = chi*rte_pll + xi*rth_pll
        ! make Raman tensor for prr process
        rtensor_pll(iel,imode,ik,iq,i,j) = rtensor_pll(iel,imode,ik,iq,i,j) + rt_pll
     END DO
  END DO
  !
  !           
  END SUBROUTINE raman_tensor
  !
END MODULE calc_raman_intensity
    
  
