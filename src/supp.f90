!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE el_phon2
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  SAVE
  !
  LOGICAL :: elph_epc
  REAL(DP) :: kx, ky, kz
  !
END MODULE el_phon2
!
!-----------------------------------------------------------------------
SUBROUTINE elphfil_epc(iq)
  !-----------------------------------------------------------------------
  !! Original routine written by Georgy Samsonidze.
  !! Rewritten by Nguyen 08/10/2022 to write el-ph matrix elements
  !! which is subsequently processed by the QERaman code.
  !-----------------------------------------------------------------------
  USE cell_base, ONLY : ibrav, alat, omega, tpiba, at, bg
  USE disp, ONLY : nq1, nq2, nq3, nqs, x_q, wq, lgamma_iq
  USE dynmat, ONLY : dyn, w2
  USE el_phon, ONLY : el_ph_mat, el_ph_mat_rec, &
       el_ph_mat_rec_col, done_elph
  USE fft_base, ONLY : dfftp, dffts, dfftb
  USE gvect, ONLY : ngm_g, ecutrho
  USE io_global, ONLY : ionode, ionode_id
  USE ions_base, ONLY : nat, nsp, atm, ityp, tau
  USE kinds, ONLY : DP
  USE klist, ONLY : xk, wk, nelec, nks, nkstot, ngk
  USE lsda_mod, ONLY : nspin, isk
  USE modes, ONLY : nirr, nmodes, npert, npertx, u, t, tmq, &
       name_rap_mode, num_rap_mode
  USE lr_symm_base, ONLY : irgq, nsymq, irotmq, rtau, gi, gimq, &
       minus_q, invsymq
  USE mp, ONLY : mp_bcast, mp_sum
  USE mp_images, ONLY : intra_image_comm
  USE mp_pools, ONLY : npool, intra_pool_comm
  USE qpoint, ONLY : xq, nksq, nksqtot, ikks, ikqs, eigqts
  USE start_k, ONLY : nk1, nk2, nk3, k1, k2, k3
  USE symm_base, ONLY : s, invs, ft, nrot, nsym, nsym_ns, &
       nsym_na, ft, sr, sname, t_rev, irt, time_reversal, &
       invsym, nofrac, allfrac, nosym, nosym_evc, no_t_rev
  USE wvfct, ONLY : nbnd, et, wg
  USE gvecw, ONLY : ecutwfc
  USE io_files, ONLY : prefix
  USE constants, ONLY : rytoev, degspin

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iq

  INTEGER :: iuelph, ios, irr, ii, jj, kk, ll, ikk, ikq, ik, &
       nu, mu, iii, n, ipert, jpert, ibnd, jbnd, vu
  CHARACTER :: cdate*9, ctime*9, sdate*32, stime*32, &
       stitle*32, myaccess*10, mystatus*7
  CHARACTER(LEN=80) :: filelph

  REAL(DP), PARAMETER :: ryd2mev  = rytoev * 1.0E3_DP
  REAL(DP), PARAMETER :: eps = 0.01/ryd2mev
  REAL(DP) :: w, weight, wspin, factor
  REAL(DP), ALLOCATABLE :: xk_collect(:,:), wk_collect(:)
  REAL(DP), ALLOCATABLE :: et_collect(:,:), wg_collect(:,:)
  INTEGER, ALLOCATABLE :: ngk_collect(:)
  INTEGER, ALLOCATABLE :: ikks_collect(:), ikqs_collect(:)
  COMPLEX(DP), ALLOCATABLE :: el_ph_mat_collect(:,:,:,:)
  !
  COMPLEX(DP) :: el_ph_sum(3 * nat, 3 * nat)
  COMPLEX(DP) :: el_ph_sum_aux(3 * nat, 3 * nat)
  REAL(DP) :: gamma, g2, w_1, w_2
  COMPLEX(DP) :: gamma2
  REAL(DP) :: epc(nbnd, nbnd, nksqtot, 3 * nat)
  COMPLEX(DP) :: epc2(nbnd, nbnd, nksqtot, 3 * nat)
  !
  INTEGER :: ftau(3,48)
  INTEGER, EXTERNAL :: find_free_unit, atomic_number

  filelph = TRIM(prefix) // '.elph'

  DO irr = 1, nirr
     IF (.NOT. done_elph(irr)) RETURN
  ENDDO

  IF (iq .EQ. 1) THEN
     myaccess = 'sequential'
     mystatus = 'replace'
  ELSE
     myaccess = 'append'
     mystatus = 'old'
  ENDIF
  ! parallel case: only first node writes
  IF (ionode) THEN
     iuelph = find_free_unit()
     OPEN(unit = iuelph, file = TRIM(filelph), FORM = 'formatted', & 
                 ACCESS = myaccess, STATUS = mystatus, iostat = ios)
  ELSE
     iuelph = 0
  ENDIF
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore('elphfil_epc', 'opening file ' // filelph, ABS(ios))

  IF (iq .EQ. 1) THEN
     CALL cryst_to_cart(nqs, x_q, at, -1)
     ! write header
     IF (ionode) THEN
         WRITE(iuelph, '(5x, a)') ' Electron-phonon matrix elements M(k,q) = sqrt(hbar/2*omega)<psi(k+q,j)|dvscf_q*psi(k,i)>'
         WRITE(iuelph, '(5x, a)') ' nbnd   nmodes   nqs   nkstot   nksqtot'
         WRITE(iuelph, '(5i9)') nbnd, nmodes, nqs, nkstot, nksqtot
         WRITE(iuelph, '(5x, a)') ' qx   qy   qz   weight_q   iq'
         WRITE(iuelph, '(5x, a)') ' kx   ky   kz   weight_k   ik   ik+q'
         WRITE(iuelph, '(5x, a)') ' ibnd  jbnd  imode  enk[eV]  enk+q[eV]  omega(q)[meV]   |M|[meV]   Re(M)[meV]   Im(M)[meV]'
         WRITE(iuelph, '(5x, a)') REPEAT('-', 78)
     ENDIF
     CALL cryst_to_cart(nqs, x_q, bg, 1)
  ENDIF

  ! collect data for current q-point
  ALLOCATE(xk_collect(3, nkstot))
  ALLOCATE(wk_collect(nkstot))
  ALLOCATE(et_collect(nbnd, nkstot))
  ALLOCATE(wg_collect(nbnd, nkstot))
  ALLOCATE(ngk_collect(nkstot))
  ALLOCATE(ikks_collect(nksqtot))
  ALLOCATE(ikqs_collect(nksqtot))
  ALLOCATE(el_ph_mat_collect(nbnd, nbnd, nksqtot, nmodes))
  IF (npool > 1) THEN
     CALL poolcollect(3, nks, xk, nkstot, xk_collect)
     CALL poolcollect(1, nks, wk, nkstot, wk_collect)
     CALL poolcollect(nbnd, nks, et, nkstot, et_collect)
     CALL poolcollect(nbnd, nks, wg, nkstot, wg_collect)
     CALL ipoolcollect(1, nks, ngk, nkstot, ngk_collect)
     CALL jpoolcollect(1, nksq, ikks, nksqtot, ikks_collect)
     CALL jpoolcollect(1, nksq, ikqs, nksqtot, ikqs_collect)
     CALL el_ph_collect(nmodes, el_ph_mat, el_ph_mat_collect, nksqtot, nksq)
  ELSE
     xk_collect(1:3, 1:nks) = xk(1:3, 1:nks)
     wk_collect(1:nks) = wk(1:nks)
     et_collect(1:nbnd, 1:nks) = et(1:nbnd, 1:nks)
     wg_collect(1:nbnd, 1:nks) = wg(1:nbnd, 1:nks)
     ngk_collect(1:nks) = ngk(1:nks)
     ikks_collect(1:nksq) = ikks(1:nksq)
     ikqs_collect(1:nksq) = ikqs(1:nksq)
     el_ph_mat_collect(1:nbnd, 1:nbnd, 1:nksq, 1:nmodes) = &
          el_ph_mat(1:nbnd, 1:nbnd, 1:nksq, 1:nmodes)
  ENDIF
  CALL cryst_to_cart(nkstot, xk_collect, at, -1)
  ! write data for current q-point
  IF (ionode) THEN
     !
     !IF (nspin .eq. 1) THEN
     !  wspin = 1.0d0 / degspin
     !ELSE
     !  wspin = 1.0d0
     !ENDIF
     !
     WRITE(iuelph, '(4f12.7, i9)') x_q(:,iq), wq(iq), iq
     !
     el_ph_sum(:,:) = (0.0d0, 0.0d0)
     DO ik = 1, nksqtot
       ikk = ikks(ik)
       ikq = ikqs(ik)
       WRITE(iuelph, '(4f12.7, 2i9)') xk_collect(:,ikk), wk_collect(ikk), ikk, ikq
       !
       DO ibnd = 1, nbnd
         DO jbnd = 1, nbnd
           !see Eq. A2 arXiv:cond-mat/0504077v2 
           DO jpert = 1, 3 * nat
             DO ipert = 1, 3 * nat
               ! note: elphmat(j,i,k,mu)=<psi_{k+q,j}|dvscf_q(mu)*psi_{k,i}>
               ! el_ph_mat_collect(j,i,k,mu) = elphmat(j,i,k,mu)
               el_ph_sum(ipert,jpert) = CONJG(el_ph_mat_collect(jbnd, ibnd, ik, ipert)) &
                      * el_ph_mat_collect(jbnd, ibnd, ik, jpert) 
             ENDDO
           ENDDO
           CALL symdyn_munu_new(el_ph_sum, u, xq, s, invs, rtau, irt, &
               at, bg, nsymq, nat, irotmq, minus_q)
           !
           DO nu = 1, nmodes
             IF (w2(nu) > 0.d0) THEN
               w = SQRT(w2(nu))
               factor = 1 / (2.0d0 * SQRT(w2(nu)))
             ELSE
               w = 0.0d0
               factor = 0.0d0
             ENDIF
             gamma = 0.0d0
             DO mu = 1, nmodes
               DO vu = 1, nmodes
                 gamma = gamma + DBLE(CONJG(dyn(mu, nu)) * &
                      el_ph_sum(mu, vu) * dyn(vu, nu))
               ENDDO
             ENDDO
             gamma = DSQRT(gamma * factor)
             epc(jbnd, ibnd, ik, nu) = epc(jbnd, ibnd, ik, nu) + gamma
             !
             gamma2 = 0.d0
             DO ipert = 1, 3 * nat
                gamma2 = gamma2 + dyn(ipert,nu) * el_ph_mat_collect(jbnd, ibnd, ik, ipert)
             ENDDO
             gamma2 = gamma2 * DSQRT(factor)
             epc2(jbnd, ibnd, ik, nu) = gamma2
             !
             WRITE(iuelph, '(3i9, 2f12.4, f20.10, 3e20.10)') ibnd, jbnd, nu, &
                rytoev * et_collect(ibnd, ikk), rytoev * et_collect(jbnd, ikq), &
                ryd2mev * w, ryd2mev * epc(jbnd, ibnd, ik, nu), &
                ryd2mev * REAL(epc2(jbnd, ibnd, ik, nu)), ryd2mev * AIMAG(epc2(jbnd, ibnd, ik, nu))
             !
           ENDDO
           !
         ENDDO
       ENDDO
       !
     ENDDO
     !
     CLOSE (unit = iuelph, status = 'keep')
     !
  ENDIF
  CALL cryst_to_cart(nkstot, xk_collect, bg, 1)
  DEALLOCATE(xk_collect)
  DEALLOCATE(wk_collect)
  DEALLOCATE(et_collect)
  DEALLOCATE(wg_collect)
  DEALLOCATE(ngk_collect)
  DEALLOCATE(ikks_collect)
  DEALLOCATE(ikqs_collect)
  DEALLOCATE(el_ph_mat_collect)
  !
  RETURN
  !
END SUBROUTINE elphfil_epc

