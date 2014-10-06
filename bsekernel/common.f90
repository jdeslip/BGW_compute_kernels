module common_m

  implicit none

  public ! only types in this module

!> Maximum number of bands supported by the *inread* routines. This sets the
!! size of arrays such as "occupations". These arrays should all be allocated
!! dynamically in the future.
  integer, parameter :: MAX_BANDS = 1000000 ! "occupations" array => 7MB
!> Maximum number of {k,q}-points supported by the *inread* routines.
!! The actual number of k-points/q-points in the WFN/bsemat/epsmat files
!! can be larger.
  integer, parameter :: MAX_KPTS = 100000 ! "kpt_read" array => 0.8 MB

!> parameters for real-space resolution in cell-truncation schemes
  integer, parameter :: n_in_box = 2
  integer, parameter :: n_in_wire = 4

!> parameter for construction of Wigner-Seitz cell
  integer, parameter :: ncell = 3

!> number of Monte-Carlo integration points
  integer, parameter :: nmc_coarse = 250000
  integer, parameter :: nmc_fine = 2500000
  integer, parameter :: nmc = nmc_fine

!> type definitions following the convention of Numerical Recipes
!! do not ever use single-precision!!
!  integer, parameter :: SP = kind(1.0)
  integer, parameter :: DP = kind(1.0d0)
!  integer, parameter :: SPC = kind((1.0,1.0))
  integer, parameter :: DPC = kind((1.0d0,1.0d0))

!> a shift on the grid in order to avoid the singularity for truncation
  real(DP), parameter :: trunc_shift(3) = (/0.5d0, 0.5d0, 0.5d0/)

!> physical constants
!!
!! These are the "2010 CODATA recommended values" taken from
!! "The NIST Reference on Constants, Units, and Uncertainty"
!! http://physics.nist.gov/cuu/
!!
!! The following variables are used throughout the package:
!!     'BOHR', 'bohr' is Bohr radius, in Angstrom
!!     'RYD', 'ryd2eV', 'rydberg' is Rydberg constant times hc, in eV
!!     'HARTREE', 'hartree' is Hartree energy, in eV
!!     'LIGHTSPEED' is inverse alpha (fine-structure constant)
!!
!! These variables are defined in the following files:
!!     Common/nrtype.f90
!!     Common/svninfo.f90
!!     Common/wfn_utils.cpp
!!     MeanField/EPM/ff2vq.py
!!     MeanField/EPM/sysParams.f90
!!     MeanField/EPM/vca.py
!!     MeanField/ICM/icm.cpp
!!     Visual/common.py
!!
  real(DP), parameter :: BOHR = 0.52917721092_dp
  real(DP), parameter :: RYD = 13.60569253_dp
  real(DP), parameter :: LIGHTSPEED = 137.035999074_dp

!> mathematical constants
!!  real(SP), parameter :: PI_S = 3.1415926535897932384626433832795_sp
  real(DP), parameter :: PI_D = 3.1415926535897932384626433832795_dp
  real(DP), parameter :: TOL_Small = 1.0d-6
  real(DP), parameter :: TOL_Zero = 1.0d-12
  real(DP), parameter :: TOL_Degeneracy = 1.0d-6
  real(DP), parameter :: INF = 1.0d12

!> string constants appearing at output
  character(len=16) :: ondate = 'Run on'
  character(len=16) :: attime = 'at'
!---------------------------

  type crystal
    real(DP) :: celvol !< cell volume in real space (a.u.)
    real(DP) :: recvol !< cell volume in reciprocal space (a.u.)
    real(DP) :: alat !< lattice constant in real space (a.u.)
    real(DP) :: blat !< lattice constant in reciprocal space (a.u.)
    real(DP) :: avec(3,3) !< lattice vectors in real space (alat)
    real(DP) :: bvec(3,3) !< lattice vectors in reciprocal space (blat)
    real(DP) :: adot(3,3) !< metric tensor in real space (a.u.)
    real(DP) :: bdot(3,3) !< metric tensor in reciprocal space (a.u.)
    integer :: nat !< number of atoms
    integer, pointer :: atyp(:) !< atomic species, atyp(1:nat)
    real(DP), pointer :: apos(:,:) !< atomic positions, apos(1:3,1:nat) (alat)
  end type crystal

!---------------------------

  type kpoints
    integer :: nspinor = 1 !< nspinor = 2 if doing two-component spinor calculation; 1 is default
    integer :: nspin   !< nspin = 1 or 2; nspin = 1 when npsinor = 2
    integer :: nrk     !< number of k-points
    integer :: mnband  !< max number of bands
    integer :: nvband  !< number of valence bands
    integer :: ncband  !< number of conduction bands
    integer  :: kgrid(3) !< Monkhorst-Pack number of k-points in each direction
    real(DP) :: shift(3) !< Monkhorst-Pack shift of grid
    real(DP) :: ecutwfc            !< wave-function cutoff, in Ry
    integer, pointer :: ngk(:)     !< number of g-vectors for each k-point
    integer :: ngkmax              !< max(ngk(:))
    integer, pointer :: ifmin(:,:) !< lowest occupied band (kpoint,spin)
    integer, pointer :: ifmax(:,:) !< highest occupied band (kpoint,spin)
    real(DP), pointer :: w(:)      !< weights (kpoint) (between 0 and 1)
    real(DP), pointer :: rk(:,:)   !< k-vector (3, kpoint) in crystal coords
    real(DP), pointer :: el(:,:,:) !< band energies (band, kpoint, spin)
    real(DP), pointer :: elda(:,:,:) !< band energies before eqp correction
    real(DP), pointer :: occ(:,:,:)  !< occupations (between 0 and 1)
    integer, pointer :: degeneracy(:,:,:) !< size of deg. subspace for (band, kpoint, spin)
  end type kpoints
  
!---------------------------
  
  type symmetry
    integer :: ntran         !< number of operations in full group
    integer :: ntranq        !< number of operations in small group of q
    real(DP) :: rq(3)        !< The q-point this ntranq belongs to
    integer :: mtrx(3,3,48)  !< symmetry matrix
    real(DP) :: tnp(3,48)    !< fractional translations
    integer :: indsub(48)    !< symmetry operations in subgroup of q
    integer :: kgzero(3,48)  !< Umklapp vectors for subgroup symmetry operations
    integer :: cell_symmetry !< 0 = cubic, 1 = hexagonal
  end type symmetry
      
!---------------------------

  type grid
    integer :: nr  !< number in reduced zone
    integer :: nf  !< number in full zone
    real(DP) :: sz !< radius of a spherical subzone equivalent to
                   !! one point in the set f
    integer, pointer :: itran(:) !< sym op to go from irrbz to fullbz
    integer, pointer :: indr(:)  !< irrbz k/q-point mapped to fullbz
    integer, pointer :: kg0(:,:) !< Umklapp vectors (for Wigner-Seitz cell)
    real(DP), pointer :: r(:,:)  !< k/q-points in reduced zone
    real(DP), pointer :: f(:,:)  !< k/q-points in full zone
  end type grid

!-----------------------------------

  type gspace
    integer :: ng       !< number of G-vectors
    integer :: nFFTgridpts !< number in FFT grid = product(FFTgrid(1:3))
    real(DP) :: ecutrho !< charge-density cutoff, in Ry
    integer, pointer :: components(:,:) !< the G-vectors, in units of 2 pi / a
    integer :: FFTgrid(3)  !< gsm: FFTgrid is the size of the FFT grid, not the maximum G-vector   
    integer, pointer :: index_vec(:) ! mapping to FFT grid
    real(DP), pointer :: ekin(:) !< kinetic energy of G-vectors
  end type gspace

!---------------------------
!> Parameters for scissors operators
!> e_cor = e_in + es + edel * (e_in - e0)
!> es and e0 are in eV. edel is a dimensionless slope.
  type sub_scissors_t
    real(DP) :: es
    real(DP) :: edel
    real(DP) :: e0
  end type sub_scissors_t

  type scissors_t
    type(sub_scissors_t) :: val
    type(sub_scissors_t) :: cond
  end type scissors_t

!---------------------------

  type wavefunction
    integer :: ng
    integer :: nband
    integer :: nspin
    integer :: nspinor = 1 !< nspinor = 2 if doing two-component spinor calculation; 1 is default
    integer, pointer :: isort(:)
    complex(DPC), pointer :: cg(:,:,:)
  end type wavefunction

!---------------------------

!> For Epsilon: this is the wavefunction before unfolding the irr. BZ.
!! For BSE: ??
  type int_wavefunction
    integer :: nspin
    integer :: nspinor = 1 !< nspinor = 2 if doing two-component spinor calculation; 1 is default
    integer, pointer ::  ng(:)     !< (nk)
    integer, pointer :: isort(:,:) !< (ngmax, nk)
    integer, pointer :: cbi(:)
    !> I think this can be decommisioned if we use kp%rk instead
    real(DP), pointer :: qk(:,:)
    complex(DPC), pointer :: cg(:,:,:)
    complex(DPC), pointer :: cgk(:,:,:,:)
  end type int_wavefunction

!------------------------------------

  !> FHJ: valence WFNs for 1 particular kpt and 1 band
  !! It stores all bands in the case of real-space WFNs
  type valence_wfns
    integer :: nband  !< This is indeed the number of valence bands
    integer :: ngv    !< Number of G-vectors
    integer :: idx_kp !< Idx of current kpt in kp/kpq structure
    integer, pointer :: isort(:)
    !> (nband+ncrit,spin). Note: the nband+ncrit index is actually useless!
    real(DP), pointer :: ev(:,:)
    complex(DPC), pointer :: zv(:,:)   !< (ngv,spin)
    !> real-space wavefunction for all "local" val. bands (fft1,fft2,fft3,band)
    complex(DPC), pointer :: wfn_fft(:,:,:,:)
  end type valence_wfns

!-------------------------------

  !> FHJ: conduction WFNs for 1 particular kpt and all bands (!) the processor owns
  type conduction_wfns
    integer :: nband  !< This is actually the number of valence+conduction bands!
    integer :: ngc    !< Number of G-vectors
    integer :: idx_kp !< Idx of current kpt in kp structure
    integer, pointer :: isort(:)
    integer, pointer :: band_index(:)
    real(DP), pointer :: ec(:,:) !< (nband,nspin)
    complex(DPC), pointer :: zc(:,:)   !< (ngc*ncownactual,spin)
    !> real-space wavefunction for all "local" cond. bands (fft1,fft2,fft3,band)
    complex(DPC), pointer :: wfn_fft(:,:,:,:)
  end type conduction_wfns

!----------------------------

  !> splines knots and coefficients
  type spline_tck
    integer :: n              !< number of knots
    integer :: k              !< degree of spline (1=linear, etc.)
    real(DP), pointer :: t(:) !< position of knots
    real(DP), pointer :: c(:) !< splines coefficient
  end type spline_tck

!----------------------------

  type siginfo
    integer :: freq_dep    !< frequency dependence of the inverse dielectric matrix
    integer :: nFreq       !< number of frequencies used in full frequency calculation
    real(DP) :: dDeltaFreq !< frequency increment (eV) for polarizability energy denominator
    real(DP) :: dBrdning   !< Lorentzian broadening (eV) for polarizability energy denominator
    real(DP), pointer :: dFreqGrid(:)   !< Grid of Frequencies for Full Frequency
    complex(DPC), pointer :: dFreqBrd(:) !< Corresponding Broadenings for Full Frequency

    integer :: nSFreq    !< number of frequencies used in spectral functions
    real(DP) :: dDeltaSFreq !< frequency increment (eV) for spectral functions
    real(DP), pointer :: dSFreqGrid(:)   !< Grid of Frequencies for spectral functions

    integer :: exact_ch    !< compute the exact static CH
    integer :: fullConvLog !< logging CH convergence
    integer :: iwritecoul  !< flag to write vcoul
    real(DP) :: tol        !< energy tolerance for degeneracy
    integer :: ggpsum
    logical :: use_hdf5     !< with -DHDF5, whether or not we actually use hdf5
    logical :: use_xdat     !< use saved exchange matrix elements from file x.dat
    logical :: use_vxcdat   !< use saved exchange-correlation matrix elements from file vxc.dat
    integer :: nkn          !< number of k-points on which to calculate Sigma (from sigma.inp)
    integer :: nfreqeval
    integer :: nq0, nq1, nq !< Number of q->0 points, q/=0 points, and total number of q-points
    logical :: subsample !< whether we perform a subsampling of the voronoi cell of q=0
    real(DP), pointer :: subweights(:) !< (nq0) weight for each subsampling q-point
    integer :: nvband
    integer :: ntband
    integer :: igamma       !< nonzero if Gamma is the only q-point, 0 otherwise
    integer :: nspin
    integer :: spin_index(2)
    integer :: icutv               !< icutv encodes presence and type of truncation
    integer :: iwriteint           !< = 0 for comm_disk, = 1 for comm_mpi
    integer :: iuseremainder
    integer :: qgrid(3)
    integer :: iscreen !< what type of screening is present. 0 = semiconductor, 1 = graphene, 2 = metal
    integer :: fdf !< finite difference form for numerical derivative of Sigma
    integer, pointer :: indkn(:) !< mapping of k-points from sigma.inp to those in kp%rk from WFN files
    integer, pointer :: diag(:) !< energy bands for which Sigma diagonal matrix elements are calculated
    integer :: ndiag    !< number of bands contained in diag(:)
    integer :: noffdiag !< offdiag
    integer :: loff
    integer :: toff
    integer :: bmin
    integer :: bmax
    integer, pointer :: off1(:) !< offdiag <bra|
    integer, pointer :: off2(:) !< offdiag |ket>
    integer, pointer :: off3(:) !< offdiag energy at which to evaluate
    integer, pointer :: offmap(:,:)  !< sig%off1(ii) = sig%diag(sig%offmap(ii,1))
    integer, pointer :: band_index(:)
    real(DP) :: dw       !< finite difference spacing for numerical derivative of Sigma in eV
    real(DP) :: ecutb    !< energy cutoff of bare coulomb interaction in Ry
    real(DP) :: ecuts    !< energy cutoff of screened coulomb interaction in Ry
    real(DP) :: xfrac    !< fraction of bare exchange
    real(DP) :: gamma    !< GPP broadening
    real(DP) :: sexcut   !< GPP SX cutoff
    real(DP) :: q0vec(3)
    real(DP) :: freqevalmin
    real(DP) :: freqevalstep
    logical :: eqp_corrections !< are we using eqp.dat
    logical :: eqp_outer_corrections !< are we using eqp_outer.dat
    type(scissors_t) :: scis
    type(scissors_t) :: scis_outer
    logical  :: wfn_outer_present
    type(spline_tck) :: spl_tck       !< scissors b-spline coefficients
    type(spline_tck) :: spl_tck_outer !< scissors b-spline coeffs for WFN_outer
    real(DP) :: avgpot
    real(DP) :: avgpot_outer
    real(DP) :: truncval(3)    !< in Bohr (au)
    real(DP) :: avgcut         !< Cut in which we do cell averages on
    real(DP), pointer :: kpt(:,:)
    real(DP), pointer :: qpt(:,:) !< (3,nq) q-points in eps0mat+epsmat files, or provided in sigma.inp
    integer :: ncrit           !< number of partially occupied bands
    real(DP) :: efermi         !< Fermi level
    real(DP) :: efermi_input   !< what is set in the input file
    logical :: rfermi          !< Fermi level relative to neutral system
    complex(DPC), pointer :: vxc(:,:) !< Vxc(G)
    complex(DPC) :: wcoul0
    logical :: freplacebz
    logical :: fwritebz
    logical :: degeneracy_check_override
    logical :: offdiagsym
    logical :: qgridsym
    logical :: die_outside_sphere
    logical :: averagew
  end type siginfo

!---------------------------
  !> dielectric matrix info using comm_mpi
  type epsmpiinfo
    integer :: ngpown
    integer :: ngpown_max
    integer :: ngpown_rem
    integer, pointer :: isrtq(:,:)   !< These 3 arrays have a dimension of (1:gvec%ng,1:(nq+1)) where
    integer, pointer :: isrtqi(:,:)  !! (nq+1) is the total # of q`s including q0.
    
    integer, pointer :: igp_owner(:)
    integer, pointer :: igp_index(:)
    integer, pointer :: inv_igp_index(:)
    integer, pointer :: nmtx(:)      !< dimension of eps(q) for each q
    real(DP), pointer :: qk(:,:)
    complex(DPC), pointer :: eps(:,:,:)    !< eps(1:gvec%ng,1:ngpown,1:(nq+1))
    
    !> dimension of epsR and epsA (1:sig%nfreq,1:gvec%ng,1:ngpown,1:(nq+1))
    complex(DPC), pointer :: epsR(:,:,:,:)
    complex(DPC), pointer :: epsA(:,:,:,:)
  end type epsmpiinfo

!---------------------------

  type wfnkqmpiinfo
    integer, pointer ::  nkptotal(:)
    integer, pointer :: isort(:,:)
    integer, pointer :: band_index(:,:)
    real(DP), pointer :: qk(:,:)
    real(DP), pointer :: el(:,:,:)
    complex(DPC), pointer :: cg(:,:,:,:)
  end type wfnkqmpiinfo

!---------------------------

  type wfnkmpiinfo
    integer, pointer ::  nkptotal(:)
    integer, pointer :: isort(:,:)
    real(DP), pointer :: qk(:,:)
    real(DP), pointer :: el(:,:,:)
    real(DP), pointer :: elda(:,:,:)
    complex(DPC), pointer :: cg(:,:,:)
  end type wfnkmpiinfo

!---------------------------

  type wpgen
    real(DP) :: wpsq(2) !< square of free el plasma freq for each spin
    real(DP) :: nelec(2) !< number of electrons for each spin per cell
    complex(DPC), pointer :: rho(:,:) !< density, (ig, ispin)
  end type wpgen

!---------------------------

  type polarizability
    integer :: freq_dep !< frequency dependence of the inverse dielectric matrix
    integer :: freq_dep_method !< full frequency calculation. 0: Adler-Wiser; 1: Shishkin and Kresse 2006
    integer :: nFreq    !< number of frequencies used in full frequency calculation
    real(DP) :: dInitFreq  !< initial frequency (eV) for polarizability energy denominator
    real(DP) :: dDeltaFreq !< frequency increment (eV) for polarizability energy denominator
    real(DP) :: dBrdning   !< Lorentzian broadening (eV) for polarizability energy denominator
    real(DP), pointer :: dFreqGrid(:) !< Grid of Frequencies for Full Frequency 
    real(DP) :: dFreqStepIncrease
    real(DP) :: dFreqCutoff1
    real(DP) :: dFreqCutoff2

    integer :: nSFreq    !< number of frequencies used in spectral function
    real(DP) :: dInitSFreq  !< initial frequency (eV) for polarizability spectral function
    real(DP) :: dDeltaSFreq !< frequency increment (eV) for polarizability spectral function
    real(DP), pointer :: dSFreqGrid(:) !< Grid of Frequencies for spectral function

    real(DP) :: dSFreqStepIncrease
    real(DP) :: dSFreqCutoff1
    real(DP) :: dSFreqCutoff2

    type(scissors_t) :: scis
    logical :: eqp_corrections !< are we using eqp.dat and eqp_q.dat files
    complex(DPC), pointer :: dFreqBrd(:)  !< Corresponding Broadenings for Full Frequency
    integer :: fullConvLog !< logging pol matrix head & tail convergence
    integer :: iwritecoul !< flag to write vcoul
    integer :: nmtx
    integer, pointer :: nmtx_of_q(:)
    integer :: nq0, nq1, nq !< Number of q->0 points, q/=0 points, and total number of q-points
    logical :: subsample !< whether we have more than one q0 point (used in subsampled calculation)
    integer :: gcomm
    logical :: min_fftgrid   !< use the smallest possible fftbox
    ! FHJ: These flags control some experimental optimizations
    integer :: os_opt_ffts       !< optimizes calculation/reuse of FFTs (real-space WFNs)
    integer :: os_para_freqs     !< num. of frequencies to calculate in parallel
    integer :: os_nfreq_para  !< num. of epsilon frequencies held by any processor 
    integer :: os_nsfreq_para !< num. of spectral frequencies held by any processor 
    logical :: os_hdf5           !< use parallel IO?
    logical :: restart        !< are we restarting the calculation? Only ok with HDF5
    integer :: stop_after_qpt !< pretent the calculation was prematurely killed after this qpt (-1=don`t kill)
    !
    integer :: WFN_FFTgrid(3)!< max. size FFTgrid that holds all WFNs
    integer :: FFTgrid(3)    !< FFT grid to use (RHO or economical one)
    !!
    integer :: iwriteint            !< = 0 for comm_disk, = 1 for comm_mpi
    logical :: skip_epsilon
    logical :: skip_chi
    logical :: use_hdf5      !< with -DHDF5, whether or not we actually use hdf5
    logical :: need_WFNq     !< will we need the WFNq file? (nq0>0.and.valueq0==1.and.iqexactlyzero==0)
    integer :: iqexactlyzero !< 1 if the q->0 point is *exactly* zero and will be read from WFN; 0 otherwise
    integer :: valueq0       !< 1=semiconductor (read from WFNq); 2=metal (read from WFN)
    integer, pointer :: irow(:)
    integer, pointer :: isrtx(:)
    integer, pointer :: isrtxi(:)
    integer :: icutv               !< icutv encodes presence and type of truncation
    real(DP) :: truncval(3)   !< in Bohr (au)
    real(DP), pointer :: eden(:,:,:)
    real(DP), pointer :: qpt(:,:)
    !> FHJ: gme = <c,k|e^(-i(q+G).r)|v,k+q>, and the indices are:
    !! (nmtx, ncownactual, nvownactual, nspin, nrk, os_para_freqs)
    complex(DPC), pointer :: gme(:,:,:,:,:,:)
    complex(DPC), pointer :: chi(:,:,:)
    integer :: ncrit
    real(DP) :: efermi
    real(DP) :: efermi_input
    logical :: rfermi
    real(DP) :: ecuts    !< energy cutoff of screened coulomb interaction in Ry
    real(DP) :: ecutsExtra
!> Reference regarding retarded/advanced functions: Catalin`s thesis, Eq. (1.44)
    complex(DPC), pointer :: chiRDyn(:,:,:,:) !< Retarded polarizability
    complex(DPC), pointer :: chiADyn(:,:,:,:) !< Advanced polarizability
    complex(DPC), pointer :: chiTDyn(:,:,:,:) !< Spectral function of polarizability

    real(DP), pointer :: edenDyn(:,:,:,:,:) !< Dynamic energy denominator
    logical :: freplacebz
    logical :: fwritebz
    logical :: degeneracy_check_override
    real(DP) :: lin_denominator !< energy threshold below which to activate lin_denominator
    type(cvpair_info), pointer :: lin_edenDyn(:,:,:,:,:) !< energies and
    ! velocities for calculating linearized denominator in dynamic case
  end type polarizability


!--------------------------------

  type cvpair_info
    real(DP) :: vc(2) !< conduction band velocity
    real(DP) :: vv(2) !< valence band velocity
    real(DP) :: ec !< conduction band energy
    real(DP) :: ev !< valence band energy
    integer :: idx_kp !< kpoint index
    logical :: vltc !< ev - ec < TOL_Degeneracy
  end type cvpair_info

!--------------------------------

  type wfnkstates
    integer :: nkpt
    integer :: ndv
    integer, pointer :: isrtk(:)
    real(DP) :: k(3)
    real(DP), pointer :: ek(:,:)
    real(DP), pointer :: elda(:,:)
    complex(DPC), pointer :: zk(:,:)
  end type wfnkstates

!---------------------------

  type wfnkqstates
    integer :: nkpt
    integer, pointer :: isrtkq(:)
    real(DP), pointer :: ekq(:,:)
    complex(DPC), pointer :: zkq(:,:)
  end type wfnkqstates

!---------------------------------

!> Used in haydock/diag only (see epsdiag.f90)
  type epsinfo
    integer :: nq !< number of q-vectors stored
    real(DP) :: emax !< maximum length of the stored q-vectors
    real(DP), pointer :: q(:,:) !< (3, nq) coordinates of q-vectors
    real(DP), pointer :: eps(:) !< (nq) head of dielectric matrix at each q-vector
    complex(DPC) :: epshead !< head of dielectric matrix at q->0
    real(DP) :: q0vec(3) !< coordinates of the q->0 vector
  end type epsinfo

!------------------------------------

!> Used in haydock/diag only
!! Note that the bands in eqpv/eqpc are indexed with respect to the Fermi
!! level, i.e., eqpv(1,:,:) is the VBM, eqpv(2,:,:) is VMB-1, etc.
  type eqpinfo
    type(scissors_t) :: scis
    type(spline_tck) :: spl_tck       !< scissors spline coefficients
    real(DP), pointer :: evqp(:,:,:)
    real(DP), pointer :: ecqp(:,:,:)
    real(DP), pointer :: evqp_co(:,:,:)
    real(DP), pointer :: ecqp_co(:,:,:)
    real(DP), pointer :: evlda(:,:,:)
    real(DP), pointer :: eclda(:,:,:)
    real(DP), pointer :: evlda_co(:,:,:)
    real(DP), pointer :: eclda_co(:,:,:)
    real(DP), pointer :: evshift(:,:,:)
    real(DP), pointer :: ecshift(:,:,:)
    real(DP), pointer :: evshift_co(:,:,:)
    real(DP), pointer :: ecshift_co(:,:,:)
  end type eqpinfo

!------------------------------------

!> moments for Haydock
  type mmtsinfo
    integer :: nmax
    integer :: nmaxp
    real(DP) :: norm
    real(DP) :: vol
    real(DP), pointer :: an(:)
    real(DP), pointer :: bn(:)
  end type mmtsinfo
  
!-------------------------------------

  type xctinfo
    logical :: is_absorption         !< whether we are running the absorption code
    logical :: inteqp                !< whether we are interpolating
    logical :: is_periodic(3)        !< which dimensions are periodic
    integer :: idimensions           !< how many total periodic dimensions
    integer :: nkpt_co               !< number of kpts in the coarse grid
    integer :: nkptq_co              !< number of kpts in the q-shifted coarse grid
    integer :: nvb_co                !< number of valence bands in the coarse grid
    integer :: ncb_co                !< number of conduction bands in the coarse grid
    integer :: n1b_co                !< nvb_co for TDA calculations, nvb_co + ncb_co for non-TDA
    integer :: n2b_co                !< ncb_co for TDA calculations, nvb_co + ncb_co for non-TDA
    integer :: nspin
    integer :: qflag                 !< int
    logical :: read_kpoints
    integer :: ipar
    integer :: iscreen !< what type of screening is present. 0 = semiconductor, 1 = graphene, 2 = metal
    logical :: renorm_transf         !< renormalize the dcc/dvv interpolation transformation?
    !> Calculate kernel blocks other than (vc)->(v'c') transitions? This will 
    !! even include transitions such as (e,e)->(e',e'). In principle, we should
    !! always include these blocks if we are interpolating later on, but they
    !! are typically not important for semiconductors within TDA.
    logical :: extended_kernel
    !> If true, we extend the co/fi transf. to WFNs of different characters:
    !! |v> = \sum_n` d_vn`|n`>, where |n`> can be |v`> or |c`>
    !! If false, we restrict the character of the expansion coefficients:
    !! |v> = \sum_v` d_vv`|v`>
    logical :: unrestricted_transf
    !> Zero out dvc/dcv coefficients
    logical :: zero_unrestricted_contrib
    logical :: patched_sampling      !< simplest case of non-uniform sampling. See absorption.inp.
    logical :: patched_sampling_co   !< Use non-uniform sampling for coarse grid for Kernel. See absorption.inp.
    integer :: zero_q0_element       !< Zero q=0 matrix element of BSE Hamiltonian? See absorption.inp
    logical :: tda                   !< use Tamm-Dancoff approximation? (Absorption only)
    integer :: iabsorp0              !< 1 means noeh_only, 0 otherwise
    integer :: iwriteint             !< = 0 for comm_disk, = 1 for comm_mpi
    logical :: eqp_corrections       !< do we use eqp.dat and eqp_q.dat
    logical :: eqp_co_corrections    !< do we use eqp_co.dat
    
!> For Coulomb interaction truncation
    integer :: iwritecoul
    integer :: icutv              !< icutv encodes presence and type of truncation
    real(DP) :: truncval(3)       !< in Bohr (au)
    integer :: nint               !< number of intervals used in
                                  !! double integral truncated_factor
    logical :: use_hdf5       !< with -DHDF5, whether or not we actually use hdf5
    logical :: bLowComm       !< If this is true, each processor will store the entire epsilon matrix
    logical :: delaunay_interp!< use Delaunay interpolation?
    integer :: neps           !< Number of G vectors to capture the dielectric cutoff
    integer :: ilowmem
    logical :: skipinterp
    integer :: ivpar, icpar
    integer :: nn             !< PlotXct restrict_kpoints
    integer :: ng
    integer :: nktotal        !< total number of unit cells
    !> Number of vertices in co k-grid that are used to expand each k-point in
    !! the fine grid for the **kernel** interpolation. This is 1 for the
    !! greedy interpolation (previous behaviour of the code), and ndims+1
    !! if we are performing Delaunay interpolation.
    integer :: npts_intp_kernel
    real(DP) :: eta           !< energy resolution
    real(DP) :: sigma         !< (used to calculate the optical spectrum)
    real(DP) :: gamma         !< (used to calculate the optical spectrum)
    real(DP) :: qshift
    real(DP) :: shift(3)      !< shift vector (this is the small shift,
                              !< used to generate WFNq_fi, referenced only if xct%read_kpoints)
    real(DP) :: lpol          !< norm of pol
    real(DP) :: pol(3)        !< light polarization for transition matrix elements
    integer :: nmtxmax        !< max. number of columns in epsmat or eps0mat
    integer :: maxpet         !< max. number of eps columns a PE can have
    integer :: theory         !< theory level in kernel calculation
                              !< 0 - GW-BSE, 1 - TDDFT
    integer :: qgrid(3)
    real(DP) :: q0vec(3)      ! This is a hack for passing q0vec for
                              ! TDDFT calculations (never used otherwise)
                              ! when there is no epsilon
    real(DP) :: short_range_frac_fock  !< Short range exchange fraction
    real(DP) :: long_range_frac_fock   !< Long range exchange fraction
    real(DP) :: screening_length       !< Screening length
                                       !< The above 3 parameters are used 
                                       !< only for TDDFT calculations.
    !> (nmtxmax): for all PEs, a global eps column igp gets mapped the local
    !! column igp_l = epsown(igp)
    integer, pointer :: epsown(:)
    !> (maxpet, npes): local eps column igp_l owned by processor ipe
    !! corresponds to the global column index ig = epsowni(igp_l, ipe)
    integer, pointer :: epsowni(:,:)
    !> (npes): number of eps cols each processor owns
    integer, pointer :: maxpe(:)
    integer, pointer :: indexq(:), isrtqi(:,:), nmtxa(:)
    integer, pointer :: ifmax(:,:), ifmaxq(:,:)
    real(DP) :: ecute         !< energy cutoff used in dielectric matrix
    real(DP) :: scaling       !< multiply kernel by arbitrary factor
    real(DP) :: ecutg         !< energy cutoff used in wavefunctions and interaction
                              !< kernel, see Rohlfing & Louie, PRB 62(8),p. 4938
                              !< (must be slightly longer than xct%ecute because of umklapp vectors)
    real(DP) :: efermi        !< computed efermi
    real(DP) :: efermi_input  !< as set in input file
    logical :: rfermi         !< relative or absolute Fermi level
    complex(DPC), pointer :: epsdiag(:,:) !< (nmtxmax, nq+1)
    !> Regular comm: (nmtxmax, maxpet, nq+1). The local processor stores row ig
    !! and a "local column" igp_l from epsilon in epscol(ig, igp_l, ik).
    !! Low comm: (nmtxmax, nmtxmax, nq+1). Each processor stores all eps(0)mat.
    !! local column = global epsilon column.
    complex(DPC), pointer :: epscol(:,:,:)
    type (wpgen) :: wpg !< charge density/fxc for TDDFT

!> Used in haydock/diag only
    integer :: nkpt_fi         !< number of kpts in the fine grid
    integer :: nkptq_fi        !< number of kpts in the q-shifted fine grid
    integer :: nvb_fi          !< number of valence bands in the fine grid
    integer :: ncb_fi          !< number of conduction bands in the fine grid
    real(DP) :: avgcut         !< Cut in which we do cell averages on
    real(DP) :: wplasmon
    complex(DPC) :: wcoul0
    integer :: vmin,vmax
    integer :: rgrid(3) !< regular grid used to calculate qpt_averages (default is kgrid)
    logical :: freplacebz
    logical :: fwritebz
    logical :: degeneracy_check_override
    logical :: die_outside_sphere
    logical :: averagew
  end type xctinfo

!------------------------------------

  type flags

!>
!> Used in haydock, diag, nonlinearoptics
!>
!> Defined flags:
!>
!>  bz0  = 0 --> use symmetries to unfold the Brillouin zone in WFN_fi file
!>         1 --> do not unfold the BZ in WFN_fi file (default)
!>  bzq  = 0 --> use symmetries to unfold the BZ in WFNq_fi file
!>         1 --> do not unfold the BZ in WFNq_fi file (default)
!>  bzc  = 0 --> use symmetries to unfold the BZ in WFN_co file
!>         1 --> do not unfold the BZ in WFN_co file (default)
!>
!>  read_dtmat = false --> calculate dcc,dvv matrices (default)
!>               true  --> take dcc,dvv matrices from file dtmat
!>
!>  eig = 0  --> do not write eigenvectors (default)
!>      < 0 --> write all eigenvectors
!>      > 0 --> write the first flag%eig eigenvectors
!>
!>  read_epsdiag = false --> read files 'eps0mat'/'epsmat' (default)
!>                 true  --> read file 'epsdiag.dat'
!>
!>  krnl = 0 --> spin triplet kernel, direct kernel only (only allow for nspin = 1)
!>         1 --> spin singlet kernel (default)
!>         2 --> local-fields + RPA, exchange kernel only
!>
!>  opr = 0 --> use velocity operator
!>        1 --> use momentum operator
!>        2 --> use JDOS operator (Haydock only)
!>
!>  lor = 0 --> use Lorentzian broadening
!>        1 --> use Gaussian broadening
!>        2 --> use Voigt broadening
!>
!>  spec = 0 --> go through the whole exciton calculation (default)
!>         1 --> calculate only absorption spectrum (this option skips
!>               all calculation and goes right to the end of the code)
!>
!>  vm = 0 --> calculate velocity/momentum matrix elements (default)
!>       1 --> read velocity/momentum matrix elements from file vmtxel
!>       2 --> use vectors from previous iteration (Haydock only!)
!>
!>  job = 0 --> ultrafast calculation
!>  job = 1 --> two-photon calculation
!>

    integer :: bz0
    integer :: lor
    integer :: bzq
    integer :: bzc
    logical :: read_dtmat
    integer :: eig
    logical :: read_epsdiag
    integer :: krnl
    integer :: opr
    integer :: spec
    integer :: vm
    integer :: job
  end type flags

!---------------------------------

  type otherinfo
    integer :: ithreeD
    integer :: knx
    integer :: kny
    integer :: knz
    real(DP) :: keta
  end type otherinfo

!---------------------------------

  type windowinfo
    real(DP) :: evalue
    real(DP) :: emax
    real(DP), pointer :: cstates(:)
    real(DP), pointer :: estates(:)
    integer, pointer :: istates(:)
    integer :: nstates
  end type windowinfo

!-----------------------------
  !> coarse-grid wavefunctions for diag/haydock
  type tdistgwf
    
    integer :: ngm !< maximum number of G-vectors
    integer :: ngl !< local number of G-vectors
    integer :: tgl !< local to global translation

!> local to global index translation : ig_g = ig_l + tgl
!> ig_g = 1 ... ng(ik) is the global G-index
!> ig_l = 1 ... ngl is the local G-index

    integer :: nk !< number of k-points
    integer :: ns !< number of spin components
    integer :: nspinor = 1 !< nspinor = 2 if doing two-component spinor calculation; 1 is default
    integer :: nv !< number of valence bands
    integer :: nc !< number of conduction bands

    integer, pointer :: ng(:) !< (nk)
    integer, pointer :: isort(:,:) !< (ngl,nk)
    complex(DPC), pointer :: zv(:,:,:,:) !< (ngl,nv,ns*nspinor,nk) 
    complex(DPC), pointer :: zc(:,:,:,:) !< (ngl,nc,ns*nspinor,nk)

  end type tdistgwf

!-----------------------------
!> MJ: work arrays - getting rid of save statements

  type work_genwf
    integer :: ikold = 0
    integer :: ikoldq = 0
    integer :: nb
    integer :: ng
    integer :: ns
    integer :: nspinor = 1 !< nspinor = 2 if doing two-component spinor calculation; 1 is default
    complex(DPC), pointer :: cg(:,:,:)
    complex(DPC), pointer :: ph(:)
    integer, pointer :: ind(:)
    integer, pointer :: isort(:)
  end type work_genwf

!-----------------------------
!> (gsm) work arrays - getting rid of save statements

  type twork_scell
    integer :: dNfft(3)
    integer :: Nplane
    integer :: Nrod
    complex(DPC), pointer :: fftbox_1D(:,:)
  end type twork_scell

  !> FHJ: mean-field header
  type mf_header_t
    integer :: version
    character(len=3) :: sheader
    character(len=32) :: sdate
    character(len=32) :: stime
    integer :: iflavor
    type(crystal) :: crys
    type(kpoints) :: kp
    type(symmetry) :: syms
    type(gspace):: gvec
  end type mf_header_t

  !> FHJ: header information for kernel files (bsedmat, bsexmat, bsemat.h5)
  type kernel_header_t

    ! Mean-field and other general information 
    type(mf_header_t) :: mf !< mf header containing number of k-points, WFN cutoff, etc.
    integer :: iscreen !< screening flag
    integer :: icutv   !< truncation flag
    real(DP) :: ecute  !< epsilon cutoff

    ! Variables specific to kernel files
    integer :: nk      !< number of k-points 
    real(DP), pointer :: kpt(:,:)
    integer :: ns      !< number of spins
    integer :: nspinor = 1 !< nspinor = 2 if doing two-component spinor calculation; 1 is default
    integer :: nvb     !< number of valence bands in the coarse grid
    integer :: ncb     !< number of conduction bands in the coarse grid
    integer :: n1b     !< nvb_co if kernel_sz==1; nvb_co + ncb_co if kernel_sz=4
    integer :: n2b     !< ncb_co if kernel_sz==1; nvb_co + ncb_co if kernel_sz=4
    integer :: theory  !< 0 for GW-BSE, 1 for TD-HF, 2 for TD-DFT
    integer :: nmat    !< number of matrices in the file (1 for bsexmat, 3 for bsedmat)
    integer :: storage !< 0 if storing full matrix (only option supported now)
    !> How many transitions blocks are there in the kernel matrix?
    !! 1 for restricted TDA kernel: vc -> v`c`
    !! 2 for restricted non-TDA kernel: {vc,cv} -> {v`c`,c`v`}  [not implemented]
    !! 4 for extended kernel: {n1,n2} -> {n1`,n2`}
    integer :: nblocks

  end type kernel_header_t

contains

  integer function bse_index(ik, ic, iv, is, xct, ncband, nvband)
    integer, intent(in) :: ik, ic, iv, is
    type(xctinfo), intent(in) :: xct
    integer, optional, intent(in) :: ncband !< default is xct%ncb_fi
    integer, optional, intent(in) :: nvband !< default is xct%nvb_fi

    integer :: ncband_, nvband_

    ! The optionals are needed for the parallelization scheme sometimes, to be set to 1.
    if(present(ncband)) then
      ncband_ = ncband
    else
      ncband_ = xct%ncb_fi
    endif

    if(present(nvband)) then
      nvband_ = nvband
    else
      nvband_ = xct%nvb_fi
    endif

    bse_index = is + (iv - 1 + (ic - 1 + (ik - 1)*ncband_)*nvband_)*xct%nspin
    return
  end function bse_index

end module common_m
