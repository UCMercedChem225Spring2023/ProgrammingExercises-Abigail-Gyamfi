Program Vector_Vector
! This code computes the dot-product of two randomly generated
! vectors using the intrinsic 'dot_product function' with dimension 
! values of 1000, 2000, 4000, 8000, and 16000. The dot product is 
! a scalar (number) hence it is declared as Real. 
! Real: specifies that the matrices will store real numbers.
!
! The program also declares two real arrays, Vector1 and Vector2 of 
! dimension (:) NDim, which are used to store the two random vectors. 
! Keep in mind that for vector, the dimension (:) is used to declare an array. 
! Dimension (:,:) is used for matrix dimension (N * M). I changed the F10.4 to F10.8 
! (8 decimal places) to get 0.00000300 s not 0.0000 s
!
! Abigail Gyamfi, CHEM 225 Spring 2023
!
! Variable Declarations
!
Implicit None
Integer::NDim
Real::Time_Start,Time_End,Scalar
Real, Dimension(:), Allocatable :: Vector1, Vector2
Character(Len=64)::Character_Temp
!
! Format Statements.
!
1000 Format(1x,'Dimension=',I6,' Time: ',F10.8,' s.')
!
! Get the dimensionality of the problem from the command line.
!
Call Get_Command_Argument(1,Character_Temp)
Read(Character_Temp,*) NDim
!
! Allocate memory for Vector1 and Vector2. Then fill the vectors with
! random values.
!
Allocate(Vector1(NDim),Vector2(Ndim))
Call random_number(Vector1)
Call random_number(Vector2)
!
! Compute the dot product of 2 vectors (multiplication of two vec1 and vec2)
! using the intrinsic dot_product function.
!
Call cpu_time(Time_Start)
Scalar = dot_product(Vector1,Vector2)
Call cpu_time(Time_End)
Write(*,1000) NDim,Time_End-Time_Start
!
End Program Vector_Vector
