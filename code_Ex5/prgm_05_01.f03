      program prgm_05_01
!
!     This program carries out calculation of 1D particle-in-a-box values.
!     Specifically, the user provides the particle mass, box length, and quantum
!     numbers for two canonical PIB eigenstates.
!
!     This program is written in atomic units.
!     
!
!
!     Variable Declarations, variable b added, t(KE), t(PE)
!
      implicit none
      integer::i,j,k,NCmdLineArgs,NBasis,IError
      real::m,l,b
      Real,Dimension(:),Allocatable::Array_MatrixH,EVals,Temp_Vector
      Real,Dimension(:,:),Allocatable::MatrixH,EVecs,Temp_Matrix
      real,external::PIB_1D_Modified_Hamiltonian_Element
      integer::n1,n2
      logical::fail
      character(len=1024)::cmd_buffer
!
!     Format Statements
!
 9000 format(1x,'Expected 1 command line arguments, but found ',i2,'.')
!
       m = 1
       l = 1
       b = 1
!      b = 10 for Q.2 ! v(X) b =10
!
!
!     NB: For H (T+V): the 5 input arguments/parameters are b, m, l, n1, and n2.
!     Read in b, m, l, n1, and n2 from the command line.
!
      fail = .false.
      NCmdLineArgs = command_argument_count()
      if(NCmdLineArgs.ne.1) then
        write(*,9000) NCmdLineArgs
        fail = .true.
      endIf
      if(fail) goto 999
      call Get_Command_Argument(1,cmd_buffer)
      read(cmd_buffer,*) NBasis
!
      Allocate(Array_MatrixH((NBasis*(NBasis+1))/2))
      Allocate(MatrixH(NBasis,NBasis),Temp_Matrix(NBasis,NBasis))
      Allocate(EVals(NBasis),EVecs(NBasis,NBasis),Temp_Vector(3*NBasis))
!
!
!     Given the input parameters, evaluate the Hamiltonian between
!     particle-in-a-box eigenfunctions n1 and n2.
!
!
      k = 1 
      do j = 1, nBasis
            do i = 1, nBasis
            If(i.ge.j) then
            Array_MatrixH(k) = PIB_1D_Modified_Hamiltonian_Element(m,l,i,j,b)
            k = k+1
         endIf
            end do
      end do
!
!
      Write(*,*)' The matrix loaded (column) lower-triangle packed:'
      Call SymmetricPacked2Matrix_LowerPac(NBasis,Array_MatrixH,MatrixH)
      Call Print_Matrix_Full_Real(MatrixH,NBasis,NBasis)
      Call SSPEV('V','L',NBasis,Array_MatrixH,EVals,EVecs,NBasis,  &
        Temp_Vector,IError)
      If(IError.ne.0) then
        Write(*,*)' Failure in DSPEV.'
        STOP
      endIf
      Write(*,*)' EVals:'
      Call Print_Matrix_Full_Real(RESHAPE(EVals,(/1,NBasis/)),1,NBasis)
      Write(*,*)' EVecs:'
      Call Print_Matrix_Full_Real(EVecs,NBasis,NBasis)  
!
!
!     The end of the job...
!
  999 continue
      end program prgm_05_01

!     The function should be a Real Function named PIB_1D_Modified_Hamiltonian_Element

      real function PIB_1D_Modified_Hamiltonian_Element(m,l,n1,n2,b)
!    
!
!     This function evaluates the potential energy matrix element < n1 | H | n2 >,
!     where n1 and n2 are particle-in-a-box eigenstate labels and H is the
!     Hamiltonian operator.
!
!
!     Variable Declarations
      implicit none
      real,intent(in)::m,l,b
      integer,intent(in)::n1,n2
      real::tMatrixElement,vMatrixElement
      real,external::PIB_1D_T_Element,PIB_1D_Modified_V_Element
!
!     The case where n1=n2 is different than n1\=n2. For this reason, we use an
!     if block to separate the evaluation of the potential energy integral for
!     these two different cases.
!
      tMatrixElement = PIB_1D_T_Element(m,l,n1,n2)
      vMatrixElement = PIB_1D_Modified_V_Element(l,n1,n2,b)
      PIB_1D_Modified_Hamiltonian_Element = tMatrixElement + vMatrixElement
!
      end function PIB_1D_Modified_Hamiltonian_Element
!
!    The function should be a Real Function named PIB_1D_T_Element

      real function PIB_1D_T_Element(m,l,n1,n2)
!
!     This function evaluates the kinetic energy matrix element < n1 | T | n2 >,
!     where n1 and n2 are particle-in-a-box eigenstate labels and T is the
!     kinetic energy operator.
!
!
!     Variable Declarations
      implicit none
      real,intent(in)::l,m
      integer,intent(in)::n1,n2
      real(kind=8), parameter :: pi = 4.D0*datan(1.D0)
!
!
!     The case where n1=n2 is different than n1\=n2. For this reason, we use an
!     if block to separate the evaluation of the kinetic energy integral for
!     these two different cases.
!
      if(n1.eq.n2) then
!     NB: when n1 = n2
      PIB_1D_T_Element = (n1**2*pi**2)/(2*l**2*m)

      else
!     NB: when n1\=n2
      PIB_1D_T_Element = 0
      endIf
!
      end function PIB_1D_T_Element

!     
!     The function should be a Real Function named PIB_1D_Modified_V_Element

      real function PIB_1D_Modified_V_Element(l,n1,n2,b)
!
!     This function evaluates the potential energy matrix element < n1 | V | n2 >,
!     where n1 and n2 are particle-in-a-box eigenstate labels and V is the
!     potential energy operator.
!
!
!     Variable Declarations
      implicit none
      real,intent(in)::l,b
      integer,intent(in)::n1,n2
      real(kind=8), parameter :: pi = 4.D0*datan(1.D0)
      real::prefactor
!
!     The case where n1=n2 is different than n1\=n2. For this reason, we use an
!     if block to separate the evaluation of the potential energy integral for
!     these two different cases.
!
      if(n1.eq.n2) then
!     NB: when n1 = n2
      PIB_1D_Modified_V_Element = (b*l)/(2)
      else
!     NB: when n1\=n2
      PIB_1D_Modified_V_Element = ((b*l)/(pi**2))* & 
      ((-1+cos(pi*(n1-n2))/((n1-n2)**2)) & 
      +(pi*sin(pi*(n1-n2)) /(n1-n2)) &
      +(1-cos(pi*(n1+n2)))/((n1+n2)**2) & 
      -(pi*sin(pi*(n1+n2))/(n1+n2)))
      endIf
!
      end function PIB_1D_Modified_V_Element
!     
      Subroutine SymmetricPacked2Matrix_LowerPac(N,ArrayIn,AMatOut)
!
!     This subroutine accepts an array, ArrayIn, that is (N*(N+1))/2 long.
!     It then converts that form to the N-by-N matrix AMatOut taking
!     ArrayIn to be in lower-packed storage form. Note: The storage mode
!     also assumes the lower-packed storage is packed by columns.
!
      Implicit None
      Integer,Intent(In)::N
      Real,Dimension((N*(N+1))/2),Intent(In)::ArrayIn
      Real,Dimension(N,N),Intent(Out)::AMatOut
!
      Integer::i,j,k
!
!     Loop through the elements of AMatOut and fill them appropriately from
!     Array_MatrixH.
!
!
! *************************************************************************
! WRITE CODE HERE TO READ THE ARRAY ELEMENTS FROM THE INPUT FILE.
! *************************************************************************
!     code written below: using -lblas -llapack
!     gfortran -lblas -llapack -o prgm_03_03.exe prgm_03_03.f03
!     
      k = 1 
      do j = 1, N
            do i = j, N
            AMatOut(i,j) = ArrayIn(k)
            k = k+1
            end do
      end do

      do i = 1, N
            do j = i+1, N
            AMatOut(i,j) = AMatOut(j,i)
            end do
      end do
! *************************************************************************
!
!
      Return
      End Subroutine SymmetricPacked2Matrix_LowerPac
!
      Subroutine Print_Matrix_Full_Real(AMat,M,N)
!
!     This subroutine prints a real matrix that is fully dimension - i.e.,
!     not stored in packed form. AMat is the matrix, which is dimensioned
!     (M,N).
!
!     The output of this routine is sent to unit number 6 (set by the local
!     parameter integer IOut).
!
!
!     Variable Declarations
!
      implicit none
      integer,intent(in)::M,N
      real,dimension(M,N),intent(in)::AMat
!
!     Local variables
      integer,parameter::IOut=6,NColumns=5
      integer::i,j,IFirst,ILast
!
 1000 Format(1x,A)
 2000 Format(5x,5(7x,I7))
 2010 Format(1x,I7,5F14.6)
!
      Do IFirst = 1,N,NColumns
        ILast = Min(IFirst+NColumns-1,N)
        write(IOut,2000) (i,i=IFirst,ILast)
        Do i = 1,M
          write(IOut,2010) i,(AMat(i,j),j=IFirst,ILast)
        endDo
      endDo
!
      Return
      End Subroutine Print_Matrix_Full_Real
