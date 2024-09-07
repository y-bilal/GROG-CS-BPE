function [im]=grid_Fessler(data,dcf,G,N);

% For ESMRMB Parallel Imaging Workshop 2010 in Wuerzburg, Germany
% Written by Nicole Seiberlich on 06/23/2010
% Please contact nes30@case.edu with comments, questions, suggestions.
% Uses Fessler toolbox :  http://www.eecs.umich.edu/~fessler/code/


mask=true([N N]);
for u=1:size(data,1)
    imagetemp = G'*(dcf.*data(u,:).');
    im(u,:,:) = embed(imagetemp,mask);
end
