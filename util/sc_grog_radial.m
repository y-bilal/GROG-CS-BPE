function [imnew,kspave,wsx,wsy]=sc_grog_radial(sig,kx,ky,stepsize,os,wsx,wsy);

% This is the basic 2D GROG code; it shifts the data in sig along
% trajectory kx-ky onto a Cartesian grid
%
% Assumptions:  Sig must be multi-channel data!  These channels must
% contain sensitivity variations in both the x and y directions (i.e. at
% least 6 coils are needed), works best with more channels equidistantly
% spaced around object.
%
% For more details see:
%  Seiberlich N, Breuer FA, Blaimer M, Barkauskas K, Jakob PM, Griswold MA.
%  Non-Cartesian data reconstruction using GRAPPA operator gridding (GROG).
%  Magn Reson Med. 2007 Dec;58(6):1257-65.

%
% Try [imnew,kspave]=sc_grog_radial(sig,kx,ky);
%
% Inputs:   sig     -- Non-Cartesian signal to be regridded (coils x read
% points x projections)
%           kx      -- kx trajectory (read points x projections)
%           ky      -- ky trajectory (read points x projections)
%           stepsize-- GRAPPA Operator step size (0.1 is usually OK)
%           os      -- oversampling factor of base Cartesian data (usually 1)
% Outputs:  kspave  -- Cartesian k-space
%           imnew	 -- Cartesian image
% 
% 
%  Written by Nicole Seiberlich 03/07/06




matsize=ceil(max(max(max(kx))-min(min(kx)),max(max(ky))-min(min(ky))))+2;   % calculate size of final matrix + 1

if nargin<5
    os=1;
end

if nargin<4
    stepsize=0.1;
end

tic;

[nc,rsize,psize]=size(sig);

regr=zeros(nc,matsize*os,matsize*os);
point=zeros(nc,matsize*os,matsize*os);
kspave=zeros(nc,matsize*os,matsize*os);

% This portion of the code calculates the weights
% It can be swapped out for the calibration along other trajectories
% The rest of the code is not trajectory dependent
if nargin<7
[wsy,wsx]=SelfCalibratingWeightsRadial(kx,ky,sig);  % Self-Calibrating radial weights
end


% Pre-calculate small shift weights in kx and ky directions 
cnt=0;
[V,D]=eig(wsy);
wssy=V*(D.^(stepsize/os))/V;
for u=-.5/stepsize:1:.5/stepsize
    cnt=cnt+1;
    wssmally(cnt,:,:)=wssy^u;
end

cnt=0;
[V,D]=eig(wsx);
wssx=V*(D.^(stepsize/os))/V;
for u=-.5/stepsize:1:.5/stepsize
    cnt=cnt+1;
    wssmallx(cnt,:,:)=wssx^u;
end


for r=1:rsize
    for p=1:psize
        diffy=(os*ky(r,p)-round(os*ky(r,p)));   % Calculate y difference
        diffx=(os*kx(r,p)-round(os*kx(r,p)));   % Calculate x difference
        numstepy=(round(diffy/stepsize));       % Calculate difference in stepsize units
        numstepx=(round(diffx/stepsize));       % Calculate difference in stepsize units
        % Find new Cartesian point by applying appropriate weight set
        fixtxy=squeeze(wssmally(.5/stepsize+1+numstepy,:,:))*squeeze(wssmallx(.5/stepsize+1+numstepx,:,:))*sig(:,r,p);        
        % Add newly calculated point to any points previously placed at
        % Cartesian location
        regr(:,round(os*ky(r,p))+round((matsize-1)/2*os)+1,round(os*kx(r,p))+round((matsize-1)/2*os)+1)=(fixtxy)+regr(:,round(os*ky(r,p))+round((matsize-1)/2*os)+1,round(os*kx(r,p))+round((matsize-1)/2*os)+1);  
        % Calculate how many points have been placed in this Cartesian
        % location
        point(:,round(os*ky(r,p))+round((matsize-1)/2*os)+1,round(os*kx(r,p))+round((matsize-1)/2*os)+1)=point(:,round(os*ky(r,p))+round((matsize-1)/2*os)+1,round(os*kx(r,p))+round((matsize-1)/2*os)+1)+1;
    end
end

% Divide summed k-space by number of points shifted to each location, as
% effective DCF
kspave(point>0)=regr(point>0)./point(point>0);
kspave=kspave(:,1:end-1,1:end-1);
% kspave=kspave(:,1:end-2,1:end-2);
% Perform IFFT to get image
imnew=permute(ifftshift(ifft2(ifftshift(permute(kspave,[2 3 1])))),[3 1 2]);
% Do Sum-of-Squares to combine coil images
t=toc;disp(['Reconstruction GROG:  ' num2str(t)]);