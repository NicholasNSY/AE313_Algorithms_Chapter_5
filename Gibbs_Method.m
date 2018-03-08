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
% 1. Calculate r1, r2, and r3: r = sqrt(sum(rvec.*rvec)) = sqrt(X^2+Y^2+Z^2)
% 2. Calculate: Cvec12 = cross(rvec1,rvec2)
%               Cvec23 = cross(rvec2,rvec3)
%               Cvec31 = cross(rvec3,rvec1)
% 3. Verify that postion vectors are coplanar: Uhatr1.*Chat23 = 0
% 4. Calculate Nvec, Dvec, and Svec: 
%    Nvec = (r1.*Cvec23)+(r2.*Cvec31)+(r3.*Cvec12)
%    Dvec = Cvec12+Cvec23+Cvec31
%    Svec = rvec1.*(r2-r3)+rvec2.*(r3-r1)+rvec3.*(r1-r2)
% 5. Calculate Velocity Vector 2: 
%    vvec2 = (sqrt(mu/(N*D))).*(((cross(Dvec,rvec2))./r2)+Svec)
% 6. Use Position Vector 2 (rvec2) and Velocity Vector 2 (vvec2) to compute
%    the orbital elements by means of ALGORITHM 4.2
%    (Orbital_Elements_from_State_Vector.m) 
%    
% Note: Sequence 5 and 6 can be calculated with either one(1) of the three(3)
%       position or velocity vectors. 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% ALGORITHM_5.1
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clc
clear
close
tic

fprintf("ALGORITHM 5.1 (Gibbs_Method.m)\n\n");

rvec1x = -294.32;
rvec1y = 4265.1;
rvec1z = 5986.7;

% rvec1x = input('r vector X value: ');       % X,Y,Z values for r vector 1
% rvec1y = input('r vector Y value: ');
% rvec1z = input('r vector Z value: ');
% fprintf('\n');

rvec1 = [rvec1x,rvec1y,rvec1z];               % r1 vector matrix

rvec2x = -1365.5;
rvec2y = 3637.6;
rvec2z = 6346.8;

% % rvec2x = input('r vector X value: ');     % X,Y,Z values for r vector 2
% % rvec2y = input('r vector Y value: ');
% % rvec2z = input('r vector Z value: ');
% % fprintf('\n');

rvec2 = [rvec2x,rvec2y,rvec2z];               % r2 vector matrix

rvec3x = -2940.3;
rvec3y = 2473.7;
rvec3z = 6555.8;

% % rvec3x = input('r vector X value: ');    % X,Y,Z values for r vector 3
% % rvec3y = input('r vector Y value: ');
% % rvec3z = input('r vector Z value: ');
% % fprintf('\n');

rvec3 = [rvec3x,rvec3y,rvec3z];              % r3 vector matrix

r1 = sqrt(sum(rvec1.*rvec1));                % Magnitude of r1 vector
r2 = sqrt(sum(rvec2.*rvec2));                % Magnitude of r2 vector
r3 = sqrt(sum(rvec3.*rvec3));                % Magnitude of r3 vector

Cvec12 = cross(rvec1,rvec2);                 % Obtain Cvec
Cvec23 = cross(rvec2,rvec3);
Cvec31 = cross(rvec3,rvec1);

Uhatr1 = rvec1/sqrt(sum(rvec1.*rvec1));      % Coplanar Check Sequence
Chat23 = Cvec23/sqrt(sum(Cvec23.*Cvec23));

CC = dot(Uhatr1,Chat23);                     % CC approx. 0 = coplanar

fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
fprintf('Coplanar Check Sequence value = %.4f\n', CC);
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n');

Nvec = (r1.*Cvec23)+(r2.*Cvec31)+(r3.*Cvec12);  % N vector matrix
N = sqrt(sum(Nvec.*Nvec));                      % Magnitude of N vector

fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
display(Nvec);
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n');

Dvec = Cvec12+Cvec23+Cvec31;                     % D vector matrix
D = sqrt(sum(Dvec.*Dvec));                       % Magnitude of D vector

fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
display(Dvec);
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n');

Svec = rvec1.*(r2-r3)+rvec2.*(r3-r1)+rvec3.*(r1-r2);   % S vector matrix

fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
display(Svec);
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n');

mu = 398600;                                  % GP for Earth

vvec2 = (sqrt(mu/(N*D))).*(((cross(Dvec,rvec2))./r2)+Svec);%v2 vector matrix

fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
display(vvec2);
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n');

beep
toc                                           % End                        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% NicholasNSY (2018)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fprintf('\n')
fprintf('Kappa KappaGold KappaPride?\n')      % Kappa KappaGold KappaPride?