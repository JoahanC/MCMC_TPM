gfortran -o ../WISE-steep-MCMC-PJDFC ../WISE-steep-MCMC-PJDFC.f ../Nedlib.a
echo 23606
echo H=18.27 from mpcorb.s3m
../WISE-steep-MCMC-PJDFC << LAST
18.27,0.2,0,23606
57231.361,216.817,-23.313,0.36555,15.33,+0.09,12.50,+0.08, 9.99,+9.99, 9.99,+9.99
57340.448,324.809, +9.925,0.91317,17.06,+0.40,14.37,+0.10, 9.99,+9.99, 9.99,+9.99
LAST
