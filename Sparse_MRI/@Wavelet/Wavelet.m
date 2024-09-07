function res = Wavelet(filterType, filterSize, wavScale)
% res = Wavelet(Filtertype, filterSize, wavScale)
%
% implements a wavelet operator
%
% (c) Michael Lustig 2007
[filterType,TI]=strtok(filterType,'_');
if isempty(TI)==0
    res.TI = 1;
else
    res.TI = 0;
end
res.adjoint = 0;
res.qmf = MakeONFilter(filterType, filterSize);
res.wavScale = wavScale;
res = class(res,'Wavelet');
