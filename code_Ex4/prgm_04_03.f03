Program Matrix_VectorProduct
! This code computes the product of a matrix and a vector to yield a new vector.
! dimension values of 1000, 2000, 4000, 8000, and 16000
! Vector1 and Vector2 are declared as an array using Real, Dimension(:)
! Real: specifies that the matrix and vector will store real numbers.
! 
! Abigail Gyamfi, CHEM 225 Spring 2023
!
! Variable Declarations
!
Implicit None
Integer::NDim
Real::Time_Start,Time_End
Real, Dimension(:,:), Allocatable :: Matrix1
Real, Dimension(:), Allocatable :: Vector1, Vector2
Character(Len=64)::Character_Temp
!
! Format Statements.
!
1000 Format(1x,'Dimension=',I6,' Time: ',F10.4,' s.')
!
! Get the dimensionality of the problem from the command line.
!
Call Get_Command_Argument(1,Character_Temp)
Read(Character_Temp,*) NDim
!
! Allocate memory for Matrix1 and Vector1. Then fill each array with
! random values.
!
Allocate(Matrix1(NDim, NDim), Vector1(NDim), Vector2(NDim))
Call random_number(Matrix1)
Call random_number(Vector1)
!
! Compute the product of a matrix and a vector using the intrinsic MatMul function.
!
Call cpu_time(Time_Start)
Vector2 = MatMul(Matrix1,Vector1)
Call cpu_time(Time_End)
Write(*,1000) NDim,Time_End-Time_Start
!
End Program Matrix_VectorProduct