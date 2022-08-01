gfortran -o WISE-steep-MCMC-PJDFC WISE-steep-MCMC-PJDFC.f Nedlib.a
echo 05693
echo H=16.79 from mpcorb.s3m
../WISE-steep-MCMC-PJDFC << LAST
16.79,0.2,0,05693
57744.197,180.113, +2.738,0.22973,14.06,+0.08,11.42,+0.03, 9.99,+9.99, 9.99,+9.99
LAST
