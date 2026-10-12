function y = resample_about(x, fs, p, q, f0)
% assumes the spectrum of x lies in [f0-fs/2, f0+fs/2)
% signal must be column vector
n = (0:size(x,1)-1).';

xs = x .* exp(-1j*2*pi*f0*n/fs);      % move assumed band centre to DC

ys = resample(xs, p, q);

fs2 = fs*p/q;  
m = (0:size(ys,1)-1).';

y  = ys .* exp( 1j*2*pi*f0*m/fs2);     % put it back; abs(y) unaffected

end