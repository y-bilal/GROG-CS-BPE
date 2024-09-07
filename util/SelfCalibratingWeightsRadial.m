function [Gy,Gx,dist]=SelfCalibratingWeightsRadial(kx,ky,sig)

% This function calculates the base kx and ky weight sets for GROG using a radial
% trajectory
% For more details see:
%  Seiberlich N, Breuer F, Blaimer M, Jakob P, Griswold M.  
%  Self-calibrating GRAPPA operator gridding for radial and spiral trajectories.  
%  Magn Reson Med. 2008 Apr;59(4):930-5.
%
% Assumptions:  Read points are equidistant throughout radial trajectory
%               Read points less than 1DeltaK apart from one another (read
%               oversampling OK)
%               Can be used for resorted constant angular velocity spiral
%               data (please check first and last points from resorting!)

% Try [Gx,Gy]=SelfCalibratingWeightsRadial(kx,ky,sig);
%
% Inputs:   kx      -- kx trajectory (read points x projections)
%           ky      -- ky trajectory (read points x projections)
%           sig     -- Non-Cartesian signal to be used for calibration (coils x read
%                       points x projections)
% Outputs:  Gy      -- base weights in ky direction
%           Gx      -- base weights in kx direction
% 
% 
%  Written by Nicole Seiberlich 09/05/06


[nc,read,ang]=size(sig);

% for each projection
for a=1:ang
    dist(:,a)=[ky(round(read/2),a)-ky(round(read/2)-1,a),kx(round(read/2),a)-kx(round(read/2)-1,a)];    %Calculate the (constant) ky and kx difference between each point along the proj
    src  = squeeze(sig(:,1:end-1,a));       % Source points for calibration along proj
    targ = squeeze(sig(:,2:end,a));         % Target points for calibration along proj
    angws = src(:,:)*pinv(targ(:,:));    % Pseudo-inverse (as in GRAPPA) for calibration (with  Tikhonov regularization) 
    logangws(:,:,a)=logm(angws);            % Take matrix log of weights for separation into Gx and Gy
end

% Separation into log(Gx) and log(Gy) in a coil-by-coil fashion
for c1=1:nc
for c2=1:nc
    logG(:,c1,c2)=squeeze(logangws(c1,c2,:)).'*pinv(dist);  % Using distance, separate into entries of weight sets
end
end

% Take matrix exponential to get Gx and Gy
Gy=expm(squeeze(logG(1,:,:)));
Gx=expm(squeeze(logG(2,:,:)));

