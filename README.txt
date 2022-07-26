The thermophysical model can run for spherical models, or for tri-axial
ellipsoids [slower].

The TPM can be run fast using table interpolation or slow.  Fast is
better.  But you have to run a table generating program first. This
program is make-table-steep.f for the steep versions or make-table.f
for the fixed surface roughness versions.

The TPM could have a variable maximum slope on the craters, leading
to variable surface roughness.  But this would make fast running
impossible.  So there are versions with a fixed maximum slope (45
degrees), such as WISE-rc-MCMC, which is not fast, and versions
which use a 75 degree maximum slope with a variable fraction of
flat terrain which leads to a variable average or RMS roughness,
such as WISE-steep-MCMC-PJDFC.f.

The PJDFC stands for the parameters in addition to the rotation pole
and visual albedo:
P is the rotation period
J is the thermal inertia sqrt{kappa*rho*C}
D is the diameter
F is the flat terrain fraction
C is the color, p_IR/p_V.

The oldest versions used Theta_1 as a parameter.  Theta is the
dimensionless sqrt{kappa*rho*C}*sqrt{2\pi/P}/[F_sun/T_ss].  T_ss
is the sub-solar temperature with zero thermal inertia.  Since
T_sun and T_ss depend on the distance from the Sun, Theta1 is
defined at 1 AU from the Sun.  But our knowledge of priors applies
more to P and J, so these are used in the later versions.  The P
prior depends on D, and the J prior depends on P, so these need to
be separate parameters.  And fairly often we know P.

The versions that do shapes are PJDFCS with two new parameters, the
axis ratios b/a and c/b.  They also have new input fields for the
observed lightcurve amplitudes.

The MCMC chain is written out on fort.21 while the maximum posterior
model is written on fort.22 in a format designed for making spectrum plots.
The posterior is a combination of the priors and the goodness of
fit.  The chain file on fort.21 is read by a program read-WISE-rc-MCMC.f,
read-WISE-rc-MCMC-PJDFC.f, or read-WISE-rc-MCMC-PJDFCS.f.  These do not
care whether the fast or slow version generated the chain.

Here is an example of running a spherical model, 
WISE-steep-MCMC-PJDFC-Icarus.csh:

echo Icarus
../WISE-steep-MCMC-PJDFC << LAST > WISE-steep-MCMC-PJDFC-Icarus.out
16.9,0.2,2.2726,Icarus
55342.832,343.387,-24.253,1.6398,20.,9.,18.694,-1.58,10.747,0.089,8.773,0.383
55452.244,262.2603,-65.5381,0.6127,15.386,0.062,12.631,0.020,7.625,0.096,4.,9.
56795.899,329.075,-13.7543,0.6350,15.513,0.167,13.320,0.041,8.,9.,5.,9.
57191.021,199.7629,29.1482,0.0589,10.247,0.02,7.514,0.02,8.,9.,5.,9. 
57252.744,240.731,-29.153,1.1539,18.,9.,15.469,0.33,8.,9.,5.,9.
LAST
cp fort.21 ../Icarus-5.JDFC
cp fort.22 ../Icarus-5.22
cp fort.21 ~/Dropbox/Icarus-5.JDFC
cp fort.22 ~/Dropbox/Icarus-5.22

The first input line to the program is: H,sigma(H),period,name
If sigma(H) > 2, there is no optical data.
In this case the period is known.  Put 0.0 for the period if it is not
known.

The remaining input lines are one per input containing:
MJD,RA,DEC,Delta,W1,sig(W1),W2,sig(W2),W3,sig(W3),W4,sig(W4)
Note that sigmas greater than 8.99 correspond to missing data.
In the example above the first line is 4-band cryo, the next is
3-band cryo, while the rest are 2-band post-cryo.

The fitting is done in flux space even though the input is in magnitudes.
The measured flux can be negative, so the conventions are
input m = mzp-2.5*log_10(|F|) and
sigma_m = sgn(F)*2.5*log_10(1+sigma(F)/|F|)
where F and sigma(F) are the flux and sigma in DNs and mzp is the
magnitude for 1 DN.  This can lead to negative magnitude sigma's when 
the flux estimate is negative, as it is for W2 in the first line.

The program calculates where the Earth was as the time of observation,
and then adds the vector defined by RA,Dec,Delta to get the asteroid
position for geometry calculations.

The programs read-WISE-rc-MCMC-PJDFC and read-WISE-rc-MCMC-PJDFCS
read the chain files, and output some statistics.  The different
numbers of parameters mean the appropriate version of the "read"
program should be used.

The programs are all fortran, and compile on gfortran on a Mac.
Example:
gfortran -o WISE-steep-MCMC-PJDFCS WISE-steep-MCMC-PJDFCS.f READ_TABLE-steep.f ~/Ned.a

The programs all use routines from the "Ned" library.
The FORTRAN code for these are in Nedlib-mac.zip.

There is a parallel set of programs for NEOCAM which only have 2 bands.
Example: NEOCAM-steep-MCMC-PJDFC.f
 
The parameters are all controlled by prior distributions which are
evaluated by adding penalties to the "chi^2".  The actual priors
are:

For diameter, a log uniform between 1 m and 1000 km.  The actual
parameter used in the MCMC is log(D) since D cannot be negative,
so no penalty is needed until the diameter wanders out of this
range.  These boundary penalties are 10*ln(D/boundary)^2 in chi^2.

For the visual albedo, the prior is the mixture model of two Rayleigh
distributions from Wright etal (2017).  The actual parameter used is
ln(p_V) since p_V cannot be negative.  A limit at p_V = 1.5 is enforced
since if p_V > 2.5, the factor (1-A) goes negative.

For the rotation pole position, the prior is uniform in 4\pi.  Some
authors have claimed that the rotation pole tends to be perpendicular
to the orbit plane, due to the YORP effect, but my programs do not
assume this.

For rotation period, the prior is figured in terms of equatorial
rotation velocity, veq.  For D < 200 m, a log-Cauchy centered on
veq = 6.5 cm/sec is used, with a half-width at half maximum of 1
unit in the natural log.  Note that the tensile stress goes like
veq^2 so this corresponds to bodies held together by tensile strength,
but it is wrong to say these are strong solid boulders.  The tensile
strength needed to support veq = 6.5 cm/sec is akin to the strength
of a lump of flour that needs to broken up by sifting.  A strong solid
boulder could spin thousands of times faster.

For bodies with D > 200 m, the spin barrier at 2 hr period is imposed
by a strong penalty, 100*[ln(P/2 hr)]**2 added to chi^2 for P < 2 hr.
In addition, the central veq is adjusted to veq0 = 6.5(D/200 m)^{5/6}.
For D = 100 km, this gives a typical period of 7 to 8 hours. 
The actual parameter in the MCMC is ln(P).

For the thermal inertia, J = sqrt(kappa*rho*C), the prior is a log
uniform between JLO=2.5 and JHI=2500 MKS units.  If J is outside
the range JLO..JHI, a modest penalty of [ln(J/limit)]^2 is added
to chi^2.  If the period is less than 2 hours, centrifugal forces
will remove any loose regolith, so the lower limit on J is raised
to JLO=250.  The actual parameter used in the MCMC is ln(J).

The discontinuities in the P and J prior at 200 m lead to the D
prior not being uniform across this gap, so a correction term 
tabulated in PRE is added to the chi^2 to even out the D
distribution.

Note that the program can be run with no data to verify the priors,
and they run very fast when the thermophysical models don't have
to be calculated.

Surface roughness is parametrized by the crater fraction, fc.
The prior is uniform on 0..1.  The actual parameter used is
y such that fc = exp(y)/(1+exp(y)), or y = logit(fc), where
logit(x) = ln(x/(1-x)).  In order to make the prior uniform
in fc, a penalty of -2ln[4*fc(1-fc)] is applied to chi^2.

The infrared to visible albedo ratio, p_IR/p_V, is taken to be a
log normal, with the parameter being ln(p_IR/p_V) with prior 0.563+/-0.340.

Finally for the shape model, only a triaxial ellipsoid is considered.
The axis ratios b/a and c/b are the new parameters, both with ranges
0 to 1, but the prior is not uniform in b/a or c/b.  Instead, a
uniform prior is assumed for (b/a)^3 [prior median b/a = 0.80] and
uniform for (c/b)^4 [prior median c/b = 0.85].  Note that lightcurve
amplitudes give data on the distribution of b/a but there is very
little information about the c/b distribution.  For triaxial models,
the diameter is computed as the volume equivalent diameter
D = 2(abc)^{1/3}.  The actual parameters used in the MCMC are
y1 = logit[(b/a)^3] and y2 = logit[(c/b)^4] with penalties of
-2ln[4/([1+exp(-y1)][1+exp(y1)])] and -2ln[4/([1+exp(-y2)][1+exp(y2)])]
added to chi^2 to make the prior uniform.

All of the parameters have legal ranges of -infinity to +infinity
except for the pole position which is encoded as RA and Dec in
radians.  The program uses the affine sampler used in the MCMC
hammer but with two adaptations for the pole position.  The affine
sampler takes two "walkers" in the simplex, and then tries a new
parameter set given by pt = p1+(p2-p1)*z with z chosen in the range
[0.5..2] with sqrt(z) uniformly distributed.  This requires a legal
range -infinity to +infinity so it has to be modified for the pole
position.  So the program calculates the angle theta between pole2
and pole1, and then rotates pole1 in the plane defined by pole2 &
pole1 by an angle z*theta to give the trial pole position.  Then
the calculation of the volume element, normally done using z^{np-1}
where np is the number of parameters, is changed to
z^{np-2}*|sin(theta*z)/sin(theta)|.

The data fitting uses a "robust" chi^2.  For each datapoint,
the contribution to chi^2 is given by f[(model-obs)/sigma]
with f(x) = x^2 if |x|<2 or 4*|x|-4 if |x|>2.  This essentially
Winsorizes the data.  
