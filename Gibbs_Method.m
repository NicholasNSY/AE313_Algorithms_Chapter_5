%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ALGORITHM_5.1 (Gibbs_Method.m)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Nicholas Ngo Syuan Yaw (ERAU)
% AE313 02DB
% Credits: Prof. Howard D. Curtis (ERAU)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Given rvec1, rvec2, and rvec3. Use Gibbs method of preliminary orbit 
% determination.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Guide:
% 1. Calculate r1, r2, and r3: r = sqrt(X^2+Y^2+Z^2)
% 2. Calculate Cvec12 = cross(rvec1,rvec2)
%              Cvec23 = cross(rvec2,rvec3)
%              Cvec31 = cross(rvec3,rvec1)
% 3. Verify that Uhat(r1).*Cvec23 = 0
% 4. Calculate Nvec, Dvec, and Svec using: 
% 5. Calculate vvec2 using:
% 6. Use rvec2 and vvec2 to compute the orbital elements by means of 
%    ALGORITHM 4.2 (Orbital_Elements_from_State_Vector.m)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% ALGORITHM_5.1
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clc
clear
close
tic

fprintf("ALGORITHM 5.1 (Gibbs_Method.m)\n\n");