      Program prgm_01_03
!
!     This program reads a 3x3 matrix from a user-provided input file. After the
!     file is opened and read, it is closed and then printed.
!
!     Abigail Gyamfi, CHEM 225 Spring 2023

      implicit none
      integer,parameter::inFileUnitA=10,inFileUnitB=11
      integer::errorFlag,i
      real,dimension(3,3)::matrixInA,matrixInB,matrixProduct
      character(len=128)::fileNameA,fileNameB
!
!
!     Start by asking the user for the name of the data file.
!
      write(*,*)' What is the name of the input data file?'
      read(*,*) fileNameA
!     write(*,*)' What is the name of the input data file?'
      read(*,*) fileNameB
      write(*,*) 

!     Open the data file and read matrixInA from that file.
!
      open(unit=inFileUnitA,file=TRIM(fileNameA),status='old',  &
        iostat=errorFlag)
      if(errorFlag.ne.0) then
        write(*,*)' There was a problem opening the input file.'
        goto 999
      endIf
      do i = 1,3
        read(inFileUnitA,*) matrixInA(1,i),matrixInA(2,i),matrixInA(3,i)
      endDo
      close(inFileUnitA)
      
!      Open the data file and read matrixInB from that file.

      open(unit=inFileUnitB,file=TRIM(fileNameB),status='old', &
              iostat=errorFlag)
      if(errorFlag.ne.0) then
        write(*,*)' There was a problem opening the input file.'
        goto 999
      endIf
      do i = 1,3
        read(inFileUnitB,*) matrixInB(1,i),matrixInB(2,i),matrixInB(3,i)
      endDo
      close(inFileUnitB)

!
!     Call the subroutine PrintMatrix to print matrixInA and matrixInB (both matrices).
!
      call PrintMatrix3x3(matrixInA)
      call PrintMatrix3x3(matrixInB)
!
!    The matrix product can be defined using MatMul as the intrinsic function 
!        Followed by: Call the subroutine PrintMatrix to print matrix product

     matrixProduct = MatMul(matrixInA,matrixInB)
     call PrintMatrix3x3(matrixProduct)

 999 continue
      End Program prgm_01_03


      Subroutine PrintMatrix3x3(matrix)
!
!     This subroutine prints a 3x3 real matrix. The output is written to StdOut.
!
      implicit none
      real,dimension(3,3),intent(in)::matrix
      real,dimension(3,3) :: StdOut
      integer::i
!
!     Format statements.
!
 1000 format(3(2x,f5.1))
!
!     Do the printing job; in this code, printing of the 2 matrices.
!
      write(*,*)' Printing Matrix'
      do i = 1,3
              StdOut(3,i) = matrix(3,i)
              StdOut(2,i) = matrix(2,i)
              StdOut(1,i) = matrix(1,i)
      endDo
      write(*,*)
      write(*,1000)transpose(StdOut)
      write(*,*)
      return
      End Subroutine PrintMatrix3x3
