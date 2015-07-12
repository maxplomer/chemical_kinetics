program CFDlineq
implicit none
REAL, ALLOCATABLE, DIMENSION(:,:) :: Matrix
REAL, ALLOCATABLE, DIMENSION(:) :: RHS
REAL, ALLOCATABLE, DIMENSION(:) :: x
integer :: i,j,n
  
PRINT*,"Enter size of square matrix"
READ*,n 

ALLOCATE(Matrix(n,n)) 
ALLOCATE(RHS(n))
ALLOCATE(x(n))

PRINT*,"Enter elements of matrix Row wise" 

DO i = 1, n
READ*,(Matrix(i,j), j = 1, n)
END DO

PRINT*,"Enter elements of RHS"
READ*,(RHS(i), i = 1, n)

PRINT*,"Your Matrix"
DO i = 1, n
PRINT*,(Matrix(i,j), j = 1, n)
END DO  
PRINT*,"Your RHS"
DO i = 1, n
PRINT*,RHS(i)
END DO 

call SOLVE(Matrix, RHS, n, x)

PRINT*,"Solution:"
DO i = 1, n
PRINT*,x(i)
END DO 

DEALLOCATE(Matrix)
DEALLOCATE(RHS)
DEALLOCATE(x)          
end program CFDlineq

	subroutine SOLVE(Matrix, RHS, n, x)
	implicit none
	integer :: i,j,k,n
	real :: m, sum
	real :: Matrix(n,n)
	real :: RHS(n)
	real :: x(n)

	!create upper triangle matrix
	DO k = 1,n-1
		DO i = k+1,n
			m=Matrix(i,k)/Matrix(k,k)
			DO j= k,n
				Matrix(i,j)=Matrix(i,j)-Matrix(k,j)*m
			END DO
			RHS(i) = RHS(i) - m*RHS(k)
		END DO
	END DO
	
	!perform substitution
	x(n) = RHS(n)/Matrix(n,n)
	DO i = n-1,1,-1
		sum=0.
		DO j = i+1,n
			sum = sum + Matrix(i,j)*x(j)
		END DO
		x(i) = (RHS(i) - sum)/Matrix(i,i)
	END DO
	
	return
	end subroutine SOLVE