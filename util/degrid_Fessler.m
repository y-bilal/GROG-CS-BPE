function sig=degrid_Fessler(im,G);

% For ESMRMB Parallel Imaging Workshop 2010 in Wuerzburg, Germany
% Written by Nicole Seiberlich on 06/23/2010
% Please contact nes30@case.edu with comments, questions, suggestions.
% Uses Fessler toolbox :  http://www.eecs.umich.edu/~fessler/code/


for u=1:size(im,1)
    sig(u,:) = G*im(u,:).';
end
