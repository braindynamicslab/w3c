function [pow,freqs,winstarts,par] = tukeySpectra(X, SR, varargin)
% TUKEYSPECTRA calculate power spectra in Tukey windows given a collection
% of time series in X, organized in columns.
%   [pow,freqs,winstarts,par] = tukeySpectra(X, SR, ...)
% input:
%   X: time series, a N-by-D matrix. N is the number of time points, D is
%   the number of dimensions/channels.
%   SR: sampling rate of the time series.
% output:
%   pow: power spectra a N-by-D-by-W array. W is the number of windows. 
%   freqs: Nf-by-1 frequency vector. Nf = ceil((nfft+1)/2);
%   winstarts: is the time (s) of the first sample point in each window. A
%   W-by-1 vector.
%   par: a struct, storing the parameters used for the spectra analysis. 
%   
%{
Author: Mengsen Zhang <mengsenzhang@gmail.com> 11/7/2019 ~
%}
[Npt,Nchn]=size(X);


p=inputParser();
p.addParameter('winSize', Npt,@isnumeric) % length of Tukey window
p.addParameter('overlaps', 0,@isnumeric) % overlap between consecutive windows (0-1)
p.addParameter('preprocessing','detrend',@isstr) % what to do before FFT each window ('detrend','zeromean','none')
p.addParameter('nfft',[])% how many points in the DFT of each segment of X
p.addParameter('specresolution',[]) % (Hz) coarsest acceptible spectral resolution
p.addParameter('norm',false);% whether to normalize the spectrum or not
p.parse(varargin{:})
par = p.Results;

% -- Tukey window setup
lag=floor(par.winSize*(1-par.overlaps));
taper=par.overlaps*2;
win=tukeywin(par.winSize,taper);
winchns=repmat(win,1,Nchn);
Nwin=floor((Npt-par.winSize)/lag+1);% number of windows

% -- fft setup
% decide spectral resolution 
if ~isempty(par.nfft) && ~isempty(par.specresolution)
    error('You cannot input both parameters ''nfft'' and ''specresolution''! Only one or none of them. ')
elseif isempty(par.specresolution)% no zero-padding
    par.nfft = par.winSize;        
else % smallest power of 2 to give resolution finer than par.specresolution
    par.nfft = 2^ceil(log2(SR/par.specresolution));
end
par.specresolution = SR/par.nfft;

Nfreq=ceil((par.nfft+1)/2);
freqs=linspace(0,SR/2,Nfreq)';

% -- compute power spectra
pow=zeros(Nfreq,Nchn,Nwin);
winstarts=zeros(Nwin,1);
% -- transform EEG
for winIdx=1:Nwin
    spidx=(winIdx-1)*lag+1:(winIdx-1)*lag+par.winSize;
    winstarts(winIdx)=find(spidx,1)/SR;% time of the start of the window
    % -- preprocessing
    switch par.preprocessing
        case 'detrend'
            seg=detrend(X(spidx,:));
        case 'zeromean'
            seg=detrend(X(spidx,:),'constant');
        case 'none'
            seg=X(spidx,:);
    end
    seg=seg.*winchns;% mask window on signal
    % -- FFT for EEG
    seg_ft=fft(seg,par.nfft,1)/par.winSize;
    seg_ft=seg_ft(1:Nfreq,:);
    % -- power spectrum of EEG
    seg_pow=abs(seg_ft).^2;
    seg_pow(2:end,:)=seg_pow(2:end,:)*2;
    if par.norm% normalization
        seg_pow = seg_pow./repmat(trapz(freqs,seg_pow,1),Nfreq,1);
    end
    pow(:,:,winIdx)=seg_pow;
end

end