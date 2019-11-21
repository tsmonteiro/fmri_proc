function rotmat=convert_rots_into_rotmat(xr,yr,zr)
%function rotmat=convert_rots_into_rotmat(xr,yr,zr)
% inputs in degrees

%xr=0; yr=0; zr=0;
%yr=10;
xr=xr*pi/180; yr=yr*pi/180; zr=zr*pi/180;
% all calcs assume we rotate about x, then y, then z (x-y-z convention)
rotmat=[cos(yr)*cos(zr) -cos(xr)*sin(zr)+sin(xr)*sin(yr)*cos(zr)  sin(xr)*sin(zr)+cos(xr)*sin(yr)*cos(zr); ...
cos(yr)*sin(zr) cos(xr)*cos(zr)+sin(xr)*sin(yr)*sin(zr) -sin(xr)*cos(zr)+cos(xr)*sin(yr)*sin(zr); ...
-sin(yr) sin(xr)*cos(yr) cos(xr)*cos(yr)];

