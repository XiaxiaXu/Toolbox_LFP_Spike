function S=DeNoiseWavelet(signal)
% Xiaxia, delete noise by Wavelet transform
% Input
%      signal, column

wname='db4'; % wave that uesed
N=1;         % deep of decompose

[~,segs]=size(signal);
for seg=1:segs
    xx=signal(:,seg);
[THR,SORH,KEEPAPP] = ddencmp('den','wp',xx);
S(:,seg)=wdencmp('gbl',xx, wname,N,THR,SORH,KEEPAPP) ;
end
