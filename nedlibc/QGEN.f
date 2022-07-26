C	QGEN.FOR
C
C	rotations will be stored as unit-norm quaternions
C
C	all routines that return a rotation will be declared
C	COMPLEX*16 and will have names that start with Q
C
C	conversion routines
C	SUBROUTINE Q2M(Q,A) converts quaternion into a rotation matrix.
C	Matrix is such that V(sky) = A*V(body) so it is the transpose
C	of the CSC convention but matches the COBLIB convention
C	COMPLEX*16 QMAT(A) converts matrix rotation matrix A into 
C		unit-norm quaternion
C
	COMPLEX*16 FUNCTION QGEN(A,B,C,D)
C
C	makes a quaternion out of 4 real components
C
	REAL*4 A,B,C,D
	REAL*4 R(4)
	COMPLEX*16 Q
	EQUIVALENCE (Q,R)
	R(1) = A
	R(2) = B
	R(3) = C
	R(4) = D
	QGEN = Q
	RETURN
	END
