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
  !
END MODULE el_phon2
!
!-----------------------------------------------------------------------
SUBROUTINE elph_prt_epc(iq)
  !-----------------------------------------------------------------------
  !! Written by Nguyen 08/10/2022 to write el-ph matrix elements
  !! which is subsequently processed by the QERaman code.
  !! Update by Nguyen 31/01/2026 to fix some issues.
  !-----------------------------------------------------------------------
  USE disp,         ONLY : nqs, x_q
  USE dynmat,       ONLY : dyn, w2
  USE el_phon,      ONLY : el_ph_mat, done_elph
  USE io_global,    ONLY : stdout, ionode, ionode_id
  USE kinds,        ONLY : DP
  USE klist,        ONLY : xk, wk, nks, nkstot
  USE lsda_mod,     ONLY : nspin
  USE modes,        ONLY : nirr, nmodes, u
  USE mp,           ONLY : mp_bcast
  USE mp_images,    ONLY : intra_image_comm
  USE mp_pools,     ONLY : npool
  USE qpoint,       ONLY : xq, nksq, nksqtot, ikks, ikqs
  USE symm_base,    ONLY : s, irt, invs
  USE lr_symm_base, ONLY : rtau, nsymq, irotmq, minus_q
  USE cell_base,    ONLY : at, bg
  USE ions_base,    ONLY : nat
  USE wvfct,        ONLY : nbnd, et
  USE io_files,     ONLY : prefix
  USE constants,    ONLY : rytoev, degspin

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iq
  ! LOCAL variables
  ! I/O
  INTEGER :: iuelph, ios
  CHARACTER(LEN=10) :: myaccess
  CHARACTER(LEN=7)  :: mystatus
  CHARACTER(LEN=80) :: filelph
  ! Indices
  INTEGER :: irr, ik, ikk, ikq
  INTEGER :: ibnd, jbnd, ipert, jpert
  INTEGER :: mu, nu, vu
  ! Constants & temporary matrix elements
  REAL(DP), PARAMETER :: ryd2mev = rytoev * 1.0E3_DP
  REAL(DP)    :: w, factor, gamma
  REAL(DP)    :: epc_val
  COMPLEX(DP) :: epc2_val
  ! Data collection arrays
  REAL(DP),    ALLOCATABLE :: xk_collect(:,:), wk_collect(:)
  REAL(DP),    ALLOCATABLE :: et_collect(:,:)
  INTEGER,     ALLOCATABLE :: ikks_collect(:), ikqs_collect(:)
  COMPLEX(DP), ALLOCATABLE :: el_ph_mat_collect(:,:,:,:)
  ! Work arrays
  COMPLEX(DP) :: el_ph_sum(3*nat, 3*nat)

  INTEGER, EXTERNAL :: find_free_unit

  !-----------------------------------------------------------------------
  ! Check el-ph readiness
  !-----------------------------------------------------------------------
  DO irr = 1, nirr
     IF (.NOT. done_elph(irr)) RETURN
  END DO

  filelph = TRIM(prefix)//'.elph'

  IF (iq == 1) THEN
     myaccess = 'sequential'
     mystatus = 'replace'
  ELSE
     myaccess = 'append'
     mystatus = 'old'
  END IF

  IF (ionode) THEN
     iuelph = find_free_unit()
     OPEN(iuelph, FILE=TRIM(filelph), FORM='formatted', &
          ACCESS=myaccess, STATUS=mystatus, IOSTAT=ios)
  ELSE
     iuelph = 0
  END IF

  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore('elph_prt_epc', 'opening file '//filelph, ABS(ios))

  !-----------------------------------------------------------------------
  ! Collect k/q data
  !-----------------------------------------------------------------------
  ALLOCATE(xk_collect(3,nkstot), wk_collect(nkstot))
  ALLOCATE(et_collect(nbnd,nkstot))
  ALLOCATE(ikks_collect(nksqtot), ikqs_collect(nksqtot))
  ALLOCATE(el_ph_mat_collect(nbnd,nbnd,nksqtot,nmodes))

  IF (npool > 1) THEN
     CALL poolcollect(3,nks,xk,nkstot,xk_collect)
     CALL poolcollect(1,nks,wk,nkstot,wk_collect)
     CALL poolcollect(nbnd,nks,et,nkstot,et_collect)
     CALL jpoolcollect(1,nksq,ikks,nksqtot,ikks_collect)
     CALL jpoolcollect(1,nksq,ikqs,nksqtot,ikqs_collect)
     CALL el_ph_collect(nmodes, el_ph_mat, el_ph_mat_collect, nksqtot, nksq)
  ELSE
     xk_collect = xk
     wk_collect = wk
     et_collect = et
     ikks_collect = ikks
     ikqs_collect = ikqs
     el_ph_mat_collect = el_ph_mat
  END IF

  !-----------------------------------------------------------------------
  ! Header file elph data
  !-----------------------------------------------------------------------
  IF (iq .EQ. 1) THEN
     CALL cryst_to_cart(nqs, x_q, at, -1)
     ! write header
     IF (ionode) THEN
         WRITE(iuelph, '(5x, a)') ' Electron-phonon matrix elements Mij(k,q) = sqrt(hbar/2*omega)<psi_{k+q,j}|dvscf_q*psi_{k,i}>'
         WRITE(iuelph, '(5x, a)') ' nbnd   nmodes   nqs   nkstot   nksqtot'
         WRITE(iuelph, '(5i9)') nbnd, nmodes, nqs, nkstot, nksqtot
         WRITE(iuelph, '(5x, a)') ' kx   ky   kz   weight_k   ik   ik+q'
         WRITE(iuelph, '(5x, a)') ' ibnd  jbnd  imode  enk[eV]  enk+q[eV]  omega(q)[meV]   |M|[meV]   Re(M)[meV]   Im(M)[meV]'
         WRITE(iuelph, '(5x, a)') REPEAT('-', 78)
     ENDIF
     CALL cryst_to_cart(nqs, x_q, bg, 1)
  ENDIF

  !-----------------------------------------------------------------------
  ! Main loop to write el-ph matrix elements for the current q-point
  !-----------------------------------------------------------------------
  IF (ionode) THEN
     el_ph_sum = (0.0_DP,0.0_DP)

     DO ik = 1, nksqtot
        ikk = ikks_collect(ik)
        ikq = ikqs_collect(ik)

        WRITE(iuelph,'(4f16.8,2i9)') xk_collect(:,ikk), wk_collect(ikk), ikk, ikq

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
            !
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
              epc_val = DSQRT(gamma * factor)
              !
              epc2_val = (0.d0, 0.d0)
              DO ipert = 1, 3 * nat
                 epc2_val = epc2_val + dyn(ipert,nu) * el_ph_mat_collect(jbnd, ibnd, ik, ipert)
              ENDDO
              epc2_val = epc2_val * DSQRT(factor)
              !
              WRITE(iuelph, '(3i6, 2f16.8, f20.10, 3e20.10)') &
                ibnd, jbnd, nu, rytoev * et_collect(ibnd, ikk), rytoev * et_collect(jbnd, ikq), &
                ryd2mev * w, ryd2mev * epc_val, &
                ryd2mev * REAL(epc2_val), ryd2mev * AIMAG(epc2_val)
              !
            ENDDO
            !
          ENDDO
        ENDDO
     END DO

     CLOSE(iuelph)
  END IF
  CALL cryst_to_cart(nkstot, xk_collect, bg, 1)
  DEALLOCATE(xk_collect, wk_collect, et_collect)
  DEALLOCATE(ikks_collect, ikqs_collect, el_ph_mat_collect)

END SUBROUTINE elph_prt_epc
!
