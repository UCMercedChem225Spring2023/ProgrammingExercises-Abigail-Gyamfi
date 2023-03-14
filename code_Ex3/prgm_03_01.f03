      Program prgm_03_01
!
!     This program reads a file name from the command line, opens that
!     file, and loads a packed/linearized form of a square matrix. Then,
!     the packed matrix is expanded assuming a row-packed form and printed.
!     Finally, the packed matrix is expanded as a column-packed form and
!     printed.
!
!     The input file is expected to have the leading dimension (an integer
!     NDim) of the matrix on the first line. The next NDim*NDim lines each
!     have one real number each given.
!
!     Abigail Gyamfi, CHEM 225 Spring 2023 
!
      Implicit None
      Integer,Parameter::IIn=10
      Integer::IError,NDim,i,j
      Real,Dimension(:),Allocatable::Array_Input
      Real,Dimension(:,:),Allocatable::Matrix
      Character(Len=256)::FileName
!
!     Begin by reading the input file name from the command line. Then,
!     open the file and read the input array, Array_Input.
!
      Call Get_Command_Argument(1,FileName)
      Open(Unit=IIn,File=TRIM(FileName),Status='OLD',IOStat=IError)
      If(IError.ne.0) then
        Write(*,*)' Error opening input file.'
        STOP
      endIf
      Read(IIn,*) NDim
      Allocate(Array_Input(NDim*NDim),Matrix(NDim,NDim))
!
! *************************************************************************
! WRITE CODE HERE TO READ THE ARRAY ELEMENTS FROM THE INPUT FILE.
! *************************************************************************
!     code written below: i, j entry of a matrix = no. in row i & column j
!     This do-loop reads values from the given input file and stores them 
!     into a Rank 1- 1D array called "Array_Input" with Dimension of NDim*NDim.
!     The do loop iterates over the range of values 1 to NDim*NDim
!
      do i= 1, NDim*NDim
        read(IIn,*) Array_Input(i)
      end do
! *************************************************************************    
!
!
!     Convert Array_Input to Matrix and print the matrix.
!
      Write(*,*)' The matrix expanded according to a row-wise ', &
        'linear packed format:'
      Call Packed2Matrix_RowWise(NDim,NDim,Array_Input,Matrix)
      Call Print_Matrix_Full_Real(Matrix,NDim,NDim)
      Write(*,*)' The matrix expanded according to a column-wise ', &
        'linear packed format:'
      Call Packed2Matrix_ColumnWise(NDim,NDim,Array_Input,Matrix)
      Call Print_Matrix_Full_Real(Matrix,NDim,NDim)
!
      End Program prgm_03_01


      Subroutine Packed2Matrix_ColumnWise(M,N,ArrayIn,AMatOut)
!
!     This subroutine accepts an array, ArrayIn, that is M*N long. It then
!     takes those elements and converts them to the M-by-N matrix AMatOut
!     such that the elements in ArrayIn are interpreted as a packed
!     column-wise form of the matrix AMatOut.
!
      Implicit None
      Integer,Intent(In)::M,N
      Real,Dimension(N*M),Intent(In)::ArrayIn
      Real,Dimension(M,N),Intent(Out)::AMatOut
!
      Integer::i,j,k
!
!     Loop through the elements of AMatOut and fill them appropriately from
!     Array_Input.
!
!
! *************************************************************************
! WRITE CODE HERE TO READ THE ARRAY ELEMENTS FROM THE INPUT FILE.
! *************************************************************************
!     code written below: using nested do-loops for j and i 
!     This subroutine called Packed2Matrix ColumnWise, includes an argument 
!     list of (M,N,ArrayIn,AMatOut) where M & N are input arguments giving 
!     the sizes of the two dimensions of the full matrix, ArrayIn is an 
!     input argument giving the rank-1 array giving the pack storage form 
!     of the matrix, and AMatOut is an output argument returning the full 
!     (un-packed) matrix. The conversion of ArrayIn to ArrayOut is carried 
!     out taking the linear packed storage form to be in row format.
! 
!     This do loops iterate over the ranges of values 1 to M for the outer 
!     loop (looping over j) and 1 to N for the inner loop (looping over i). 
!     For each iteration of the inner loop, the corresponding element of 
!     ArrayIn (as indexed by k) is assigned to the corresponding element of 
!     AMatOut (as indexed by i and j). After each assignment, the value 
!     of k is incresed by 1 to continue to the next element of ArrayIn.
!
      k = 1
      do j = 1, M
        do i = 1,N
       AMatOut(i,j) = ArrayIn(k)
       k = k+1
        end do
      end do
! *************************************************************************
!
      Return
      End Subroutine Packed2Matrix_ColumnWise


      Subroutine Packed2Matrix_RowWise(M,N,ArrayIn,AMatOut)
!
!     This subroutine accepts an array, ArrayIn, that is M*N long. It then
!     takes those elements and converts them to the M-by-N matrix AMatOut
!     such that the elements in ArrayIn are interpreted as a packed
!     row-wise form of the matrix AMatOut.
!
      Implicit None
      Integer,Intent(In)::M,N
      Real,Dimension(N*M),Intent(In)::ArrayIn
      Real,Dimension(M,N),Intent(Out)::AMatOut
!
      Integer::i,j,k
!
!     Loop through the elements of AMatOut and fill them appropriately from
!     Array_Input.
!
!
! *************************************************************************
! WRITE CODE HERE TO READ THE ARRAY ELEMENTS FROM THE INPUT FILE.
! *************************************************************************
!     This do loop is like the columnwise code, except that indices(i,j) of 
!     AMatOut are swapped to the rowwise (j,i) used by default in Fortran
!
      k = 1
      do j = 1, M
        do i = 1,N
       AMatOut(j,i) = ArrayIn(k)
       k = k+1
        end do
      end do
! *************************************************************************
!
!
      Return
      End Subroutine Packed2Matrix_RowWise
!
!     The general full matrix printing subroutine print_matrix_full_real.f03 
!     in github repository is copied into the the source code
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
