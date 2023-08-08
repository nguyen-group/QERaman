!
! Copyright (C) 2006-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE write_p_avg(filp, spin_component, firstk, lastk)
  !---------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE cell_base,            ONLY : at, bg, ibrav
  USE constants,            ONLY : rytoev
  USE gvect,                ONLY : ngm, g
  USE lsda_mod,             ONLY : nspin
  USE ener,                 ONLY : ef
  USE wvfct,                ONLY : et, nbnd, npwx
  USE klist,                ONLY : xk, wk, nks, nkstot, ngk, igk_k
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE uspp,                 ONLY : nkb, vkb, okvan
  USE becmod,               ONLY : bec_type, becp, calbec, &
                                   allocate_bec_type, deallocate_bec_type
  USE noncollin_module,     ONLY : noncolin, npol
  USE ldaU,                 ONLY : lda_plus_u
  USE wavefunctions, ONLY : evc
  USE io_global,            ONLY : ionode, ionode_id, stdout
  USE mp,                   ONLY : mp_bcast, mp_sum
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp_world,             ONLY : world_comm
  USE uspp_init,            ONLY : init_us_2
  !
  IMPLICIT NONE
  !
  INTEGER :: spin_component, nks1, nks2, firstk, lastk, npw, ndim
  INTEGER :: iunout, ios, ik, ibnd, jbnd, ipol, nbnd_occ
  COMPLEX(DP), ALLOCATABLE :: ppsi(:,:), ppsi_us(:,:), matp(:,:,:)
  CHARACTER (len=256) :: filp, namefile
  !
  IF (lda_plus_u) CALL errore('write_p_avg', &
                       'write_p_avg not working with LDA+U',1)
  ALLOCATE(matp(nbnd,nbnd,3))
  CALL allocate_bec_type ( nkb, nbnd, becp)

  IF (nspin==1.or.nspin==4) THEN
     nks1=max(1,firstk)
     nks2=min(nkstot, lastk)
     IF (spin_component /= 1)  &
        CALL errore('write_p_avg','incorrect spin_component',1)
  ELSEIF (nspin==2) THEN
     IF (spin_component == 1) THEN
        nks1=max(1,firstk)
        nks2=min(nks/2,lastk)
     ELSEIF (spin_component==2) THEN
        nks1=nks/2 + max(1,firstk)
        nks2=nks/2 + min(nks/2,lastk)
     ELSE
        CALL errore('write_p_avg','incorrect spin_component',1)
     ENDIF
  ENDIF

  ios = 0
  IF ( ionode ) THEN
     iunout=58
     namefile=trim(filp)
     OPEN (unit = iunout, file = namefile, status = 'unknown', form = &
          'formatted', iostat = ios)
     REWIND (iunout)
  ENDIF

  CALL mp_bcast (ios, ionode_id, world_comm)
  IF ( ios/=0 ) CALL errore ('write_p_avg', 'Opening filband file', abs (ios) )

  DO ik = nks1, nks2
     !
     !   Compute the number of occupated bands at this k point
     !
     DO ibnd = 1, nbnd
        IF (et (ibnd, ik)<=ef) nbnd_occ = ibnd
     ENDDO
     IF (nbnd_occ==nbnd) WRITE( stdout, '(5x,/,&
             &"No empty band at point ", i4,3f10.5)') &
             ik,  (xk (ipol, ik) , ipol = 1, 3)
     !
!     ALLOCATE(ppsi(npwx*npol,nbnd_occ))
!     ALLOCATE(ppsi_us(npwx*npol,nbnd_occ))
!============= Adding by Nguyen 10/10/2022 ===============================
      ALLOCATE(ppsi(npwx*npol,nbnd))
      ALLOCATE(ppsi_us(npwx*npol,nbnd))
!===========================================================================
     !
     npw = ngk(ik)
     CALL init_us_2 (npw, igk_k(1,ik), xk(1,ik), vkb)
     !
     !   read eigenfunctions
     !
     CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)

!     CALL calbec ( npw, vkb, evc, becp, nbnd_occ )
     CALL calbec ( npw, vkb, evc, becp, nbnd)

     IF (noncolin) THEN
        ndim = npwx * npol
     ELSE
        ndim = npw
     END IF
!     DO ipol=1,3
!        CALL compute_ppsi(ppsi, ppsi_us, ik, ipol, nbnd_occ, spin_component)
        ! FIXME: use ZGEMM instead of DOT_PRODUCT
!        DO ibnd=nbnd_occ+1,nbnd
!           DO jbnd=1,nbnd_occ
!              matp(ibnd-nbnd_occ,jbnd,ipol)=  &
!                   DOT_PRODUCT( evc(1:ndim,ibnd),ppsi(1:ndim,jbnd) )
!              IF (okvan) THEN
!                 matp(ibnd-nbnd_occ,jbnd,ipol)=                  &
!                      matp(ibnd-nbnd_occ,jbnd,ipol)+             &
!                      (0.d0,0.5d0)*(et(ibnd,ik)-et(jbnd,ik))*  &
!                      DOT_PRODUCT( evc(1:ndim,ibnd),ppsi_us(1:ndim,jbnd))
!              ENDIF
!           ENDDO
!        ENDDO
!     ENDDO
!============= Adding by Nguyen 10/10/2022 ===============================
     ! Dipole vector for all band
     DO ipol=1,3
        ! ppsi contains P_c^+ p | psi_ik > for the ipol cartesian coordinate
        ! ppsi_us contains the additional term required for US PP
        ! see J. Chem. Phys. 120, 9935 (2004)
        CALL compute_ppsi(ppsi, ppsi_us, ik, ipol, nbnd, spin_component)
        ! FIXME: use ZGEMM instead of DOT_PRODUCT
        DO ibnd = 1, nbnd
           DO jbnd = 1, nbnd
              matp(ibnd,jbnd,ipol)=  &
                   DOT_PRODUCT( evc(1:ndim,jbnd),ppsi(1:ndim,ibnd) )
              IF (okvan) THEN
                 matp(ibnd,jbnd,ipol)=                  &
                      matp(ibnd,jbnd,ipol)+             &
                      (0.d0,0.5d0)*(et(jbnd,ik)-et(ibnd,ik))*  &
                      DOT_PRODUCT( evc(1:ndim,jbnd),ppsi_us(1:ndim,ibnd))
              ENDIF
           ENDDO
        ENDDO
     ENDDO
!===========================================================================
     DEALLOCATE(ppsi)
     DEALLOCATE(ppsi_us)
     CALL mp_sum(matp, intra_bgrp_comm)

     IF (ionode) THEN
!        IF (ik == nks1) &
!           WRITE (iunout, '(" &p_mat nbnd=",i4,", nks=",i4," /")') &
!                 nbnd, nks2-nks1+1
!        WRITE (iunout, '(10x,3f10.6,i7)') xk(1,ik),xk(2,ik),xk(3,ik), &
!                                          nbnd_occ
!
!        DO ipol=1,3
!           WRITE (iunout, '(i3)') ipol
!           DO ibnd=nbnd_occ+1,nbnd
!              WRITE (iunout, '(5f15.8)') &
!                    (abs(matp(ibnd-nbnd_occ,jbnd,ipol))**2, jbnd=1,nbnd_occ)
!============= Adding by Nguyen 03/10/2022 ===============================
        IF (ik == nks1) THEN
        WRITE(iunout, '(5x, a)') 'Matrix elements of dipole vector D(k) = <psi(k,j)||p|psi(k,i)>'
        WRITE(iunout, '(5x, a)') 'nbnd   ndnd_occ   nks'
        WRITE(iunout, '(5x,3i5)') nbnd, nbnd_occ, nks2-nks1+1
        WRITE(iunout, '(5x, a)') 'ikx   iky   ikz   weight of k'
        WRITE(iunout, '(5x, a)') 'ibnd   jbnd   eni[eV]   enj[eV]   &
           Re(p,x)   Im(p,x)   Re(p,y)   Im(p,y)   Re(p,z)   Im(p,z)'
        WRITE(iunout, '(5x, a)') REPEAT('-', 78)
        ENDIF
        !
        WRITE(iunout, '(5x,4f10.6)') xk(1,ik),xk(2,ik),xk(3,ik), wk(ik)
        !
        DO ibnd = 1, nbnd
           DO jbnd = 1, nbnd
              WRITE(iunout, '(2i6,2f11.4,6e20.10)') ibnd,jbnd,&
                 rytoev*et(ibnd,ik),rytoev*et(jbnd,ik),&
                 REAL(matp(ibnd,jbnd,1)),AIMAG(matp(ibnd,jbnd,1)),&
                 REAL(matp(ibnd,jbnd,2)),AIMAG(matp(ibnd,jbnd,2)),&
                 REAL(matp(ibnd,jbnd,3)),AIMAG(matp(ibnd,jbnd,3))
!===========================================================================
           ENDDO
        ENDDO
     ENDIF
  ENDDO

  IF (ionode) THEN
     CLOSE(iunout)
  ENDIF

  DEALLOCATE(matp)
  !
  RETURN
END SUBROUTINE write_p_avg
