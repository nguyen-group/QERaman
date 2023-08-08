  !
  ! Module of Raman program
  ! Started by Y. Tatsumi 2016/11/16
  ! Modified by N. T. Hung 2022/10/10
  !-------------------------------------------------------------------
  MODULE raman_mod
  !-------------------------------------------------------------------
  !
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = KIND(1.d0)
  !
  PRIVATE
  SAVE
  !
  PUBLIC :: init_parameter, init_polvec, &
       prefix, outdir, fil_dvec, fil_elph, &
       sorb, circular_pol, nonpol, plot_matele_opt, plot_matele_elph, plot_raman_k, &
       nks, nksq, nbnd, nspin, npol, nmode, nqs, nel, nrs, nbnd_occ, n_kinterp, &
       gamma, gamma_raman, rs_start, rs_end, &
       k, q, wk, wq, eigv, eq, polvec, &
       raman_k, intensity_raman, rs, rtensor, rtensor_k, &
       dvec, matele_opt, matele_elph, elaser, &
       pi, im, &
       ry2ev, ev2ry, ev2cm, ry2cm
  !
  !
  CHARACTER (len=256) :: prefix
  !! Prefix of the input/output data
  CHARACTER (len=256) :: outdir
  !! Directry of output
  CHARACTER (len=256) :: fil_dvec
  !! File name of matrix element of dipole vector (from modified bands.x)
  CHARACTER (len=256) :: fil_elph
  !! File name of electron-phonon matrix element (from modified ph.x)
  !
  !
  LOGICAL :: sorb
  !! Consideration of spin-ortbit interaction (sorb = .FALSE. in this version)
  LOGICAL :: circular_pol
  !! If true, Raman spectra for circular polarizaed incident light is calculated, too
  LOGICAL :: nonpol
  !! If true, calculate for non-polarized light
  LOGICAL :: plot_matele_opt
  !! If true, generate the file of matrix element for electron-photon matrix element
  LOGICAL :: plot_matele_elph
  !! If true, generate the file of matrix element for electron-phonon matrix element
  LOGICAL :: plot_raman_k
  !! If true, generate the file of Raman matrix element for each k point
  !
  !
  INTEGER :: nks
  !! Total number of k points in scf calculation by Quantum Espresso
  INTEGER :: nksq
  !! Total number of k+q points in scf calculation by Quantum Espresso
  INTEGER :: nbnd
  !! Total number of energy bands
  !INTEGER :: ngs
  !! Total number of plane wave for wave functions
  INTEGER :: nspin
  !! Number of spin index (1:without 2:with spin orbit interaction)
  INTEGER :: npol
  !! Number of polarization vector
  !INTEGER :: nbv
  !! Number of valence band (in epw calculation)
  INTEGER :: nbnd_occ
  !! Occupation number from scf calculation
  INTEGER :: nmode
  !! Number of phonon mode
  INTEGER :: nqs
  !! Number of q points
  INTEGER :: nel
  !! Number of laser energy mesh
  INTEGER :: nrs
  !! Number of Raman shift mesh
  INTEGER :: n_kinterp
  !! Number of iteration of interpolation
  !INTEGER :: resonant_mode
  !! "0" non-resonant, "1" first-order resonant
  !
  !REAL (DP) :: alat
  !! lattice constant (atomic unit, Bohr)
  REAL (DP) :: gamma
  !! Broadning factor for Raman formula
  REAL (DP) :: gamma_raman
  !! Broadning factor to plot the Raman intensity
  REAL (DP) rs_start, rs_end
  !! Claculation range of ramanshift (eV)
  !
  REAL (DP), ALLOCATABLE :: k(:,:)
  !! wave vector k(ik,idim)
  REAL (DP), ALLOCATABLE :: q(:,:)
  !! phonon wave vector q(iq,idim)
  !REAL (DP), ALLOCATABLE :: g(:,:)
  !! reciprocal lattice vector g(ig,idim)
  REAL (DP), ALLOCATABLE :: wk(:)
  !! weight of each k point wk(ik)
  REAL (DP), ALLOCATABLE :: wq(:)
  !! weight of each q point wq(iq)
  REAL (DP), ALLOCATABLE :: eigv(:,:)
  !! Eigen energy of each band at each k point eigv(ik,ib)
  REAL (DP), ALLOCATABLE :: eq(:,:)
  !! Eigen energy of phonon at each k point eq(imode,iq)
  !REAL (DP), ALLOCATABLE :: avec(:,:)
  !! primitive lattice vector
  !REAL (DP), ALLOCATABLE :: bvec(:,:)
  !! primitive reciplocal lattice vector
  REAL (DP), ALLOCATABLE :: intensity_raman(:,:,:,:)
  !! Raman intensity: intensity(iel,ipi,ips,irs)
  !REAL (DP), ALLOCATABLE :: intensity_raman_c(:,:,:)
  !! Raman intensity for sigma+ --> sigma- or sigmai --> sigma+: intensity(iel,4:5,irs)
  REAL (DP), ALLOCATABLE :: rs(:)
  !! Raman shift rs(irs)
  REAL (DP), ALLOCATABLE :: elaser(:)
  !! Incident laser energy: elaser(iel)
  !
  COMPLEX (DP), ALLOCATABLE :: dvec(:,:,:,:)
  !! Transition dipole moment (dipole vector) dvec(ik,ib_f,ib_i,ixyz)
  COMPLEX (DP), ALLOCATABLE :: matele_opt(:,:,:,:)
  !! Optical matrix element matele_opt(ik,ibi,ibf,ip)
  COMPLEX (DP), ALLOCATABLE :: matele_elph(:,:,:,:,:)
  !! Electron-phonon matrix element matele_elph(imode,iq,ik,ibi,ibf)
  !COMPLEX (DP), ALLOCATABLE :: coef_wfc(:,:,:,:)
  !! Coefficient of wave function coef_wfc(ik,ib,ig,ispin)
  COMPLEX (DP), ALLOCATABLE :: polvec(:,:)
  !! polarization (Jones') vector polvec(npol,3) 1:x, 2:y, 3:z, 4:LCP, 5:RCP
  COMPLEX (DP), ALLOCATABLE :: raman_k(:,:,:,:,:,:)
  !! Raman matrix elements for each k point raman_k(iel,ipi,ips,imode,ik,iq)
  COMPLEX(DP), ALLOCATABLE :: rtensor(:,:,:,:,:)
  !! Raman tensor rtensor(iel,imode,iq,i,j)
  COMPLEX(DP), ALLOCATABLE :: rtensor_k(:,:,:,:,:,:)
  !! Raman tensor for each k point rtensor_k(iel,imode,ik,iq,i,j)
  !
  ! parameter
  REAL (DP) pi
  !! circumference ratio
  COMPLEX (DP) im
  !! imaginary i
  !
  ! conversion of unit
  REAL (DP) ry2ev
  !! Rydberg --> eV
  REAL (DP) ev2ry
  !! eV --> Rydberg
  REAL (DP) ev2cm
  !! eV --> cm-1
  REAL (DP) ry2cm
  !! ry --> cm-1
  !
  !
  CONTAINS
  !
  !
  SUBROUTINE init_parameter
  !
  pi = acos(-1.0d0)
  !ry2ev = 27.2114d0
  ry2eV = 13.60568d0
  ev2ry = 1.0d0 / ry2ev
  ev2cm = 8065.0d0
  ry2cm = ry2ev * ev2cm
  !
  im = CMPLX(0.0d0,1.0d0)
  !
  END SUBROUTINE init_parameter
  !
  !
  SUBROUTINE init_polvec
  !
  IF (circular_pol .EQ. .FALSE. .AND.  nonpol .EQ. .FALSE.) THEN
     npol = 3
  ELSEIF (circular_pol .EQ. .TRUE. .AND. nonpol .EQ. .FALSE.) THEN
     npol = 5
  ELSEIF (nonpol .EQ. .TRUE.) THEN
     npol = 6
  ENDIF
  !
  ALLOCATE(polvec(npol,3))
  polvec(:,:) = 0.0d0
  !ips = 1, p=(1,0,0)
  polvec(1,1) = cmplx(1.0d0, 0.0d0)
  polvec(1,2) = cmplx(0.0d0, 0.0d0)
  polvec(1,3) = cmplx(0.0d0, 0.0d0)
  !ips = 2, p=(0,1,0)
  polvec(2,1) = cmplx(0.0d0, 0.0d0)
  !polvec(2,1) = cmplx(cos(pi/3.0d0), 0.0d0)
  polvec(2,2) = cmplx(1.0d0, 0.0d0)
  !polvec(2,2) = cmplx(sin(pi/3.0d0), 0.0d0)
  polvec(2,3) = cmplx(0.0d0, 0.0d0)
  !ips = 3, p=(0,0,1)
  polvec(3,1) = cmplx(0.0d0, 0.0d0)
  polvec(3,2) = cmplx(0.0d0, 0.0d0)
  polvec(3,3) = CMPLX(1.0d0, 0.0d0)
  IF (circular_pol .eq. .TRUE.) THEN
     !ips = 4, p=(1,i,0) Left handed circular polarizaed light
     polvec(4,1) = (1.0d0/SQRT(2.0d0))*CMPLX(1.0d0, 0.0d0)
     polvec(4,2) = (1.0d0/SQRT(2.0d0))*CMPLX(0.0d0, 1.0d0)
     polvec(4,3) = CMPLX(0.0d0, 0.0d0)
     !ips = 5, p=(1,-i,0) Right-handed circular polarized light
     polvec(5,1) = (1.0d0/SQRT(2.0d0))*CMPLX(1.0d0, 0.0d0)
     polvec(5,2) = (1.0d0/SQRT(2.0d0))*CMPLX(0.0d0, -1.0d0)
     polvec(5,3) = CMPLX(0.0d0, 0.0d0)
  ENDIF
  IF (nonpol .EQ. .TRUE.) THEN
     polvec(6,1) = (1.0d0/SQRT(3.0d0))*CMPLX(1.0d0, 0.0d0)
     polvec(6,2) = (1.0d0/SQRT(3.0d0))*CMPLX(1.0d0, 0.0d0)
     polvec(6,3) = (1.0d0/SQRT(3.0d0))*CMPLX(1.0d0, 0.0d0)
  ENDIF
  !
  END SUBROUTINE init_polvec
  !
  !
  END MODULE raman_mod
