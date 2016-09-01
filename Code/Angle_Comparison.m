%Code to determine differences in angles from OIM and ang file.

%Euler Angles from OIM
phi1a = 70.9;
PHIa = 14.5;
phi2a = 293.7;

%Euler Angles from Matlab using ang file (converts to degrees since ang
%file gives them in radians)
phi1b = (4.5867)*360/(2*pi);
PHIb = (1.64059)*360/(2*pi);
phi2b = (3.56219)*360/(2*pi);

A = euler2gmat(phi1a, PHIa, phi2a);
B = euler2gmat(phi1b, PHIb, phi2b);
[Angle, axis, deltaG] = GeneralMisoCalc(A,B,'cubic');