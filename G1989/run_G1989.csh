gfortran -o ../WISE-steep-MCMC-PJDFCS WISE-steep-MCMC-PJDFCS.f Nedlib.a
echo G1989
echo H=17.35 from mpcorb.s3m
../WISE-steep-MCMC-PJDFC << LAST
17.35,0.2,0,G1989
57059.028,226.092,-13.769,0.51782,16.00,+0.05,13.91,+0.03, 9.99,+9.99, 9.99,+9.99
57443.636, 80.817, -6.305,0.60931,16.03,+0.25,14.08,+0.10, 9.99,+9.99, 9.99,+9.99
57671.255,302.430,-46.992,0.37484,14.85,+0.11,12.65,+0.22, 9.99,+9.99, 9.99,+9.99
55243.147, 67.174,-12.134,0.79551,17.24,+0.17,15.19,+0.13, 9.12,+0.07, 7.05,+0.06
LAST
