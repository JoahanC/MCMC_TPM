	SUBROUTINE GG_BAR(T1,T2,MOI,TORQUE)
C
C	computes average gravity gradient torque over time interval T1..T2
C	MOI: moment of inertia tensor in S/C axes, units kg*m^2
C	TORQUE: units N*m, in celestial coordinates
C
	REAL*8 T1,T2,T
	REAL*4 MOI(3,3)
	REAL*4 TORQUE(3),NGG(3)
	REAL*4 A(3,3),R(3,3)
	REAL*4 XYZ(3),XYZ2000(3)
	REAL*8 T01JAN03 /1092096025.52D0/
	EPOCH = 2003.0 + ((T1+T2)/2-T01JAN03)/(365.25D0*864.D2)
	CALL PREMAT(EPOCH,R)
	DO I=1,3
	  R(I,I) = 1.+R(I,I)
	ENDDO
C
	DO I=1,3
	  TORQUE(I) = 0
	ENDDO
	DO I=1,8
	  T = T1+(T2-T1)*(I-0.5)/8
	  CALL ASPECT(T,A)
	  CALL ORBIT(T,XYZ)
	  CALL VMAT(R,XYZ,XYZ2000)
	  CALL GRAVITY_GRADIENT(XYZ2000,MOI,A,NGG)
	  DO J=1,3
	    TORQUE(J) = TORQUE(J)+NGG(J)/8
	  ENDDO
	ENDDO
	RETURN
	END
	
	SUBROUTINE GRAVITY_GRADIENT(XYZ,MOI,A,TORQUE)
C
C	computes the gravity gradient torque on a satellite
C	Inputs:
C	XYZ: position of satellite wrt Earth's center in km
C	MOI: moment of inertia tensor in S/C axes, units kg*m^2
C	A: aspect matrix, V(cel J2000) = A*V(s/c)
C	Output:
C	TORQUE: units N*m, in celestial coordinates
C
	REAL*4 XYZ(3),MOI(3,3),A(3,3),TORQUE(3)
C
C	MU is G*M(Earth) in MKS units
	REAL*4 MU
	PARAMETER (MU = 403.5033E12)
C
	REAL*4 RSC(3),IR(3),GG(3)
C
C	evaluate Eq(17.31) of Wertz, N_GG = (3*MU/R^3)*(hat-R x (I*hat-R))
	CALL VMAT(A,XYZ,RSC)
	CALL NORM(RSC)
	CALL MAV(MOI,RSC,IR)
	CALL CROSS(RSC,IR,GG)
	R = 1000.*SQRT(DOT(XYZ,XYZ))
	DO I=1,3
	  GG(I) = (3*MU/R**3)*GG(I)
	ENDDO
	CALL MAV(A,GG,TORQUE)
	RETURN
	END
