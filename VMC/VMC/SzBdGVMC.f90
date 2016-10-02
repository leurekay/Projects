module crystal
    
  type atomconfig     
    !atomic info
	real*8 r(3)			                  !absolute position of the atom within unitcell, in units of primitive lattice constants.
	real*8 E, B, U                        !atomic level, spin spliting, and hubbard U	
	logical kondo, mott      	          !atomic mottness and kondoness

    !independent bonds
    integer nbond                                                  !number of independent bonds
    real*8, allocatable, dimension (:) :: J, V                     !physical spin-exchange and Coulomb interaction along the bond (independent of orbs)
	real*8, allocatable, dimension (:, :) :: bond                  !bond(1:3, 1:nbond), independent bonds connecting to the atom
    integer, allocatable, dimension (:) :: nbatom                  !nbatom(1:nbond), neighbor atoms connected to iatom by bonds
    complex*16, allocatable :: t(:), tz(:)                         !physical spin independent/dependent hoppings

    !variational parameters 
    real*8 gU                                                    !local Gutzwiller penalty potential for doublon
	real*8 muloc, mzloc                                          !local variational mu and zeeman
	complex*16 pairloc                                           !local variational singlet pairing on the atom
    complex*16, allocatable :: hop(:), hopz(:)                   !variational hoppings on bonds
    complex*16, allocatable :: pair(:), pairz(:)                 !variational pairings on bonds
  end type atomconfig

  type unitcellconfig
    integer natom                                             !number of atoms within the unitcell
	real*8, dimension (3, 3) :: abc, dualabc                  !translation vectors in real and reciprocal spaces
	
	type (atomconfig), allocatable, dimension (:) :: atom	  !atom configs
  end type unitcellconfig

end module crystal

module MPI_place
  !use MPI

  ! MPI info
  integer ithread, nthread, ierr
  integer mpi_comm_world
 
end module MPI_place

module vmcplace

  use MPI_place
  use crystal
  
  !---------------- vmc inputs ---------------------------------------------------  
	character*12  model

    type (unitcellconfig), target :: unitcell

    integer nd, ndim            !nd = irreducible spatial dim of lattice
                                !ndim = nd + 1 (the last dim used for atoms within unitcell)
    integer na1uc, nambu        !natom per unitcell and the associated nambu dimension
    integer nuc, naa            !number of unitcells and number of all atoms in the lattice  
    integer limits(3)           ! limits = (/# uc along a, # uc along b, # atoms per uc/)
    integer, allocatable :: A2R(:, :)         !A2R(1:ndim, 1:naa), mapping from atom on lattice to corrd 
    integer mnn                               !maximal number of nearest-neighbor atoms for a given atom
    integer, allocatable :: nn(:, :)          !nn(1:mnn, 1:naa), nearest-neighbors of all atoms on lattice (zero if nn bond invalid)
	integer mbond                             !maximal number of independent bonds connecting atoms 
	integer, allocatable :: nghbr(:, :)       !neighbors of atoms on lattice along independent bonds
    logical antiperiodic(2)  !antiperiodic boundary conditions along a and b 

    integer Npair            !number of Cooper pairs

    logical mott, kondolattice            !mottness and kondoness
    integer kondosites                    !number of kondosites (automated)
    logical, allocatable :: kondo(:)      !indicate whether a site is a kondo site
    real*8  U, gU, Vzz, Vcc, Vdh          !U and Jastraw factors        
    !in principle these parameters depend on atoms involved, but is applied literally for all atoms
    !for the time being (could be upgraded to atom-dependent Jastraw later)

    complex*16, allocatable :: table(:, :, :)    !table(a, b, iR2R), where a and b indicate atoms in uc, and iR2R indicate inter-uc index

    logical bruteforce                           !switch to calculate determinants by brute force
    integer nsample, nbath, nmc1sample           !controls for monte carlo

	logical optimize                             !switch to perform automatic optimization for variables user specified
	integer nvar, nfield                         !number of var's to be optimized, nfield is the number of variational mean fields
	integer optinfo                              !if optinfo = 0, only var bounds are accessed, without further processing
	                                             !if optinfo = 1, dynamic memories are allocated for the first time
												 !if optinfo > 1, dynamic memories are not allocated repeatedly, but are updated
	real*8, target, dimension (100) :: var, vmin, vmax    !vars, min and max of vars to be optimized (max number of var's is 100)
	character*6, target :: vname(100)                     !names of vars to be optimized

  
  
  !------------------- dynamic part of vmc ----------------------------------------------
    complex*16, allocatable, dimension (:, :) :: D, G       !wavefunc matrix and its inverse
    complex*16 ratio, det                                   !determinant ratio, and determinant itself (used in bruteforce mode only)
    integer, allocatable :: occup(:, :)                     !occup(1:naa, spin),  occupation/tag of electrons over all atoms

	logical LMoptimize, Newton, Conjugate                   !linear minimization scheme with two options: newton and conjugate
    integer whichtable, LM                                  !whichtable picks the correct derivative table, LM = nvar + 1, used in conjugate scheme
    complex*16, target, allocatable :: tables(:, :, :, :), Dv(:, :, :)      !derivative tables and derivative determinants used in conjugate scheme
	complex*16, allocatable, dimension (:) :: Edlog, Avdlog, dlogs, ratios  !derivatives and ratios used in newton scheme
	complex*16, allocatable, dimension (:, :) :: Svv, Evv                   !matrices used in conjugate scheme


    integer iseed, isample                         !random number generator, and number of samples sampled
    integer nattempt, naccepted                    !number of attempts and number of accepted attempts

  !------------------- vmc output  --------------------------------------------------------
    real*8, dimension(2) :: rho_s, rhoav, delrho     !diamagnetic part of the superfluid density along x and y
    real*8 acceptance, Eav, delE,  Ekav, Epav      !acceptance ratio, average E, standard error of E (average kinetic E and potential E not used yet)
    real*8, allocatable :: sdwav(:), delsdw(:), cdwav(:), delcdw(:)
    
end module vmcplace

  
program SzBdGVMC_MPI
!==========================================================================!
!    This is the MPI version of SzBdGVMC performs VMC for superconducting  !
!    state with conserved Sz                                               !
!    Copy Right:   Qiang-Hua Wang                                          !
!    Last updated:  Dec 30 / 2015, Apr 21 / 2016                           !
!==========================================================================!

use vmcplace

implicit none
integer isel, i
character*20 explain

!this line is necessary if MPI is not used
ithread = 0;  nthread = 1

!initialize MPI
call MPI_INIT( ierr )
call MPI_COMM_SIZE( MPI_COMM_WORLD, nthread, ierr )
call MPI_COMM_RANK( MPI_COMM_WORLD, ithread, ierr )

do i = 0, nthread - 1
  call MPI_barrier( mpi_comm_world, ierr )  
  if(ithread /= i)cycle
  open(10, file = 'job.txt')
  read(10, *)explain;   read(10, *)model
  read(10, *)explain;   read(10, *)isel
  if(isel /= 6)then; close(10); cycle; end if
  read(10, *)explain; read(10, *)conjugate; close(10)
end do

!print*, 'type 1/2/3/4/5/6 to perform vpoint/vline/cobyla/bisec/gene/LM:';  read*, isel
  
select case (isel)
  case (1);  call vpoint()
  case (2);  call vline()
  case (3);  call cobyla_opt()
  case (4);  call bisec_opt()
  case (5);  call gene_opt()
  case (6);  call LMoptimizer()
  case default;  stop 'unknown job selection'
end select

!release MPI
call MPI_FINALIZE( ierr )

stop
end

