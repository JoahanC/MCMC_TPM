	SUBROUTINE SPRCNPX_SET
C
	INTEGER*2 JX,JY
	COMMON /SPRCNPX/ JX(0:1023),JY(0:1023)
	DATA JX(1023)/0/
C
	SAVE /SPRCNPX/
	write (*,*) 'SPRCNPX_SET called'
	DO KPIX=0,1023
	  JPIX=KPIX
	  IX=0
	  IY=0
	  IP=1
	  DO WHILE (JPIX.NE.0)			!Break up the pixel number
	    ID=MOD(JPIX,2)			!By even and odd bits to get
	    JPIX=JPIX/2				!IX and IY
	    IX=ID*IP+IX
	    ID=MOD(JPIX,2)
	    JPIX=JPIX/2
	    IY=ID*IP+IY
	    IP=2*IP
	  ENDDO
	  JX(KPIX)=IX
	  JY(KPIX)=IY
	ENDDO
	RETURN
	END
