	COMPLEX*16 FUNCTION QMUL(P,Q)
C	PRODUCT OF TWO QUATERNIONS P*Q
	REAL*4 P(4),Q(4),R(4)
	COMPLEX*16 QQ
	EQUIVALENCE (QQ,R)
C
C	ORDER OF COMPONENTS:	i,j,k,1
C	Reminder: 1*1=1; i*i=j*j=k*k=-1
C	1*i=i;  1*j=j;  1*k=k;
C	i*j=-j*i=k;  j*k=-k*j=i; k*i=-i*k=j
C
	R(1)=P(4)*Q(1)+Q(4)*P(1)+P(2)*Q(3)-Q(2)*P(3)
	R(2)=P(4)*Q(2)+Q(4)*P(2)+P(3)*Q(1)-Q(3)*P(1)
	R(3)=P(4)*Q(3)+Q(4)*P(3)+P(1)*Q(2)-Q(1)*P(2)
	R(4)=P(4)*Q(4)-P(1)*Q(1)-P(2)*Q(2)-P(3)*Q(3)
	QMUL = QQ
	RETURN
	END
