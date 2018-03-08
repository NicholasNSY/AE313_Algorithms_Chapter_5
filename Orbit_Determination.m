%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ALGORITHM_5.4 (Orbit_Determination.m)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Nicholas Ngo Syuan Yaw (ERAU)
% AE313 02DB
% Credits: Prof. Howard D. Curtis (ERAU)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Given the range (rho), azimuth (A), angular elevation (a) together with 
% the rates rhodot, Adot, and adot relative to an earth-based tracking
% station (for which the altitude (H), latitude (phi), and local sidereal 
% time (thetaLST) are known),calculate the state vectors rvec and vvec in
% the geocentric equatorial frame.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Guide:
% 1. Using the altitude (H), latitude (phi), and local sidereal time (thetaLST) of the site. 
%    Calculate its geocentric position vector R: Rvec = ((Re+H)*cos(phi)).*[cos(thetaLST),sin(thetaLST),0]+(Re+H).*[0,0,sin(phi)]
% 2. Calculate the topocentric declination (delta): delta = asin((cos(phi)*cos(A)*cos(a))+(sin(phi)*sin(a)))
% 3. Calculate the topocentric right ascension (alpha): alpha = thetaLST-h
%
%    #Note: !~Calculate hour angle (h): acos((cos(phi)*sin(a)-sin(phi)*cos(A)*cos(a))/cos(delta))~!
%
% 4. Calculate the direction cosine unit vector rho (urhovec): urhovec = cos(delta).*[cos(alpha),sin(alpha),0]+[0,0,sin(delta)]
% 5. Calculate the geocentric position vector (rvec): rvec = Rvec+(rho.*urhovec)
% 6. Calculate the inertial velocity (Rvecdot) of the site: Rvecdot = cross(omega,Rvec)
% 
%    #Note: !~omega is the angular velocity of Earth: [0,0,72.92*10^-6] rad/s~!
% 
% 7. Calculate the declination rate (deltadot):
%    deltadot = ((-Adot*cos(phi)*sin(A)*cos(a))+(adot*((sin(phi)*cos(a))-(cos(phi)*cos(A)*sin(A)))))/cos(delta)
% 8. Calculate the right ascension rate (alphadot):
%    alphadot = (((Adot*cos(A)*cos(a))-(adot*sin(A)*sin(a))+(deltadot*sin(A)*cos(a)*tan(delta)))/((cos(phi)*sin(a))-(sin(phi)*cos(A)*cos(a))))-omega
% 9. Calculate the direction cosine rate vector (urhovecdot):
%    urhovecdot = [(-alphadot*sin(alpha)*cos(delta))-(deltadot*cos(alpha)*sin(delta)),(alphadot*cos(alpha)*cos(delta))-(deltadot*sin(alpha)*sin(delta)),deltadot*cos(delta)]
% 10. Calculate the geocentric velocity vector (vvec): vvec = Rvecdot+(rhodot*urhovec)+(rho*urhovecdot)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% ALGORITHM_5.4
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clc
clear
close

fprintf("ALGORITHM 5.4 (Orbit Determination from Angle and Range Measurements)\n\n");

Adeg = 90;
Adotdeg = 0.1130;
adeg = 30;
adotdeg = 0.05651;
rho = 2551;
rhodot = 0;
thetaLSTdeg = 300;
H = 0;
phideg = 60;

% Adeg = input('Azimuth (deg): ');                  % Angle and Range values
% Adotdeg = input('Azimuth Rate (deg/s): ');
% adeg = input('Elevation (deg): ');
% adotdeg = input('Elevation Rate (deg/s): ');
% rho = input('Range (km): ');
% rhodot = input('Range Rate (km/s): ');
% thetaLSTdeg = input('Local Sidereal Time (deg): ');
% H = input('Sea Level (km): ');
% phideg = input('Latitude (deg): ');

tic

A = Adeg*(pi/180);                                % Data Conversion to Radians
Adot = Adotdeg*(pi/180);
a = adeg*(pi/180);
adot = adotdeg*(pi/180);
thetaLST = thetaLSTdeg*(pi/180);
phi = phideg*(pi/180);

Re = 6378;                                        % Radius of Earth

Rvec = ((Re+H)*cos(phi)).*[cos(thetaLST),sin(thetaLST),0]+(Re+H).*[0,0,sin(phi)]; 
                                                  % R vector matrix
                                                  
delta = asin((cos(phi)*cos(A)*cos(a))+(sin(phi)*sin(a))); 
                                                  % Topocentric Declination

deltadeg = delta*(180/pi);                        % delta Data Conversion                    

if (0<Adeg)||(Adeg<180)                           % Hour Angle
    h = 2*pi-acos((cos(phi)*sin(a)-sin(phi)*cos(A)*cos(a))/cos(delta));
else
    h = acos((cos(phi)*sin(a)-sin(phi)*cos(A)*cos(a))/cos(delta));
end

hdeg = h*(180/pi);                                % h Data Conversion
alpha = thetaLST-h;                               % Topocentric Right Ascension
alphadeg = alpha*(180/pi);                        % a Data Conversion

urhovec = cos(delta).*[cos(alpha),sin(alpha),0]+[0,0,sin(delta)]; 
                                                  % Direction Cosine Unit Vector rho 

rvec = Rvec+(rho.*urhovec);                       % Geocentric Position Vector

fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
display(rvec)
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n');

omegavec = [0,0,72.92*10^-6];                     % Angular Velocity of Earth (vec)
omega = 72.92*10^-6;                              % Angular Velocity of Earth (mag)

Rvecdot = cross(omegavec,Rvec);                   % Inertial Velocity

deltadot = ((-Adot*cos(phi)*sin(A)*cos(a))+(adot*((sin(phi)*cos(a))-(cos(phi)*cos(A)*sin(A)))))/cos(delta);
                                                  % Declination Rate
                                                  
alphadot = (((Adot*cos(A)*cos(a))-(adot*sin(A)*sin(a))+(deltadot*sin(A)*cos(a)*tan(delta)))/((cos(phi)*sin(a))-(sin(phi)*cos(A)*cos(a))))-omega;
                                                  % Right Ascension Rate

urhovecdot = [(-alphadot*sin(alpha)*cos(delta))-(deltadot*cos(alpha)*sin(delta)),(alphadot*cos(alpha)*cos(delta))-(deltadot*sin(alpha)*sin(delta)),deltadot*cos(delta)];
                                                  % Direction Cosine Rate Vector

vvec = Rvecdot+(rhodot*urhovec)+(rho*urhovecdot); % Geocentric Velocity Vector

fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
display(vvec)
fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n');

beep
toc                                           % End                        
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% NicholasNSY (2018)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fprintf('\n')
fprintf('Kappa KappaGold KappaPride?\n')      % Kappa KappaGold KappaPride?