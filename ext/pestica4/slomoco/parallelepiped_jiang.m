function [meanvox,meanvoxz]=parallelepiped_jiang(mocor_arranged)
%function [meanvox,meanvoxz]=parallelepiped_jiang(mocor_arranged)
% mocor_arranged are the motion parameters arranged in order x,y,z trans, x,y,z rotations
% so when using the 6DOF txt output from 3dvolreg, call this as 
% parallelepiped_jiang(motion(:,[6 7 5 3 4 2]));
% as 3dvolreg outputs the motion vector as volume, z,x,y rotation, z,x,y translation, then sq diff before, after moco
%
% outputs the mean total displacement on a parallelepiped about the size of a head
% also outputs the mean out-of-plane displacement of same (so x,y- are zeroed before the summation)


% calculating total dislpacement
% define bounding box of human head, constant for constant head sizes
% note we assume rotation about center of image
% AFNI and SPM report rotation about center, FSL reports rotation about image corner
% if using mcflirt, you have to convert to center first
lowz=-40; % in mm
hiz=52;
lowx=-64;
hix=64;
lowy=-84;
hiy=84;

tdim=size(mocor_arranged,1);
raddeg=0.01745329;

for d = 1:tdim
    
    x = mocor_arranged(d,1);
    y = mocor_arranged(d,2);
    z = mocor_arranged(d,3); 
    
    xr = mocor_arranged(d,4);
    yr = mocor_arranged(d,5);
    zr = mocor_arranged(d,6);


    cz=cos(zr*raddeg);
    sz=sin(zr*raddeg);
    cx=cos(xr*raddeg);
    sx=sin(xr*raddeg);
    cy=cos(yr*raddeg);
    sy=sin(yr*raddeg);

    npix=0;
    displ=0;
    displz=0;

    for zo=[lowz hiz]

        for yo=[lowy hiy]
            
            for xo=[lowx hix]
                
                xn=(cz*cy-sz*sx*sy)*xo + (sz*cy+cz*sx*sy)*yo - (cx*sy)*zo - x;
                yn=        (-sz*cx)*xo +          (cz*cx)*yo +    (sx)*zo - y;
                zn=(cz*sy+sz*sx*cy)*xo + (sz*sy-cz*sx*cy)*yo + (cx*cy)*zo - z;

                displ=displ+sqrt((xn-xo)*(xn-xo) +(yn-yo)*(yn-yo)+ (zn-zo)*(zn-zo));
                displz=displz+sqrt((zn-zo)*(zn-zo));

                npix = npix+1;
            end
        end
    end

    displ= displ/npix;
    displz= displz/npix;
    meanvox(d) = displ;
    meanvoxz(d) = displz;
end

