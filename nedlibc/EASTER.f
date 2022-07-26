	SUBROUTINE EASTER(YEAR,MONTH,DAY)
C	for input YEAR, returns MONTH & DAY of Easter
	INTEGER YEAR
	INTEGER MONTH
	INTEGER DAY
C	from http://aa.usno.navy.mil/faq/docs/easter.php on 21-Mar-2008
	INTEGER y,c,n,k,i,j,l,m,d
	y = YEAR
	c = y/100
	n = y-19*(y/19 )
	k = (c-17)/25
	i = c-c/4-(c-k)/3+19*n+15
	i = i-30*(i/30)
	i = i-(i/28)*(1-(i/28)*(29/(i+1))*((21-n)/11))
	j = y+y/4+i+2-c+c/4
	j = j-7*(j/7)
	l = i-j
	m = 3+(l+40)/44
	d = l+28-31*(m/4)
	MONTH = m
	DAY = d
	RETURN
	END
