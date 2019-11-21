function [xr,yr,zr]=convert_rotmat_into_rots(rotmat)
%function [xr,yr,zr]=convert_rotmat_into_rots(rotmat)
% outputs in radians

%rotmat=[1 0 0; 0 1 0; 0 0 1];
%rotmat=[0.944589     -0.321023     0.0685322; 0.320728      0.947044     0.0155578; -0.0698974    0.00728447  0.997528];
% all calcs assume we rotate about x-y-z (in AFNI's 3dvolreg, I think its z-x-y)
% afni reports the above leads to zr=-18.7253  xr=-0.4174  yr=-4.0082 in degrees
xr=atan2(rotmat(3,2),rotmat(3,3));
yr=atan2(-rotmat(3,1),sqrt(rotmat(3,2)^2+rotmat(3,3)^2));
zr=atan2(rotmat(2,1),rotmat(1,1));

%function rotmat=convert_rots_into_rotmat(xr,yr,zr)
% %function rotmat=convert_rots_into_rotmat(xr,yr,zr)
% % inputs in degrees

% xr=xr*pi/180; yr=yr*pi/180; zr=zr*pi/180;
% % all calcs assume we rotate about x, then y, then z (x-y-z convention)
% rotmat=[cos(yr)*cos(zr) -cos(xr)*sin(zr)+sin(xr)*sin(yr)*cos(zr)  sin(xr)*sin(zr)+cos(xr)*sin(yr)*cos(zr); ...
% cos(yr)*sin(zr) cos(xr)*cos(zr)+sin(xr)*sin(yr)*sin(zr) -sin(xr)*cos(zr)+cos(xr)*sin(yr)*sin(zr); ...
% -sin(yr) sin(xr)*cos(yr) cos(xr)*cos(yr)];

