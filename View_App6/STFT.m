function S = STFT(x,window_length)

dx = length(x);
dw = window_length;
S = zeros(dx,dw);

if size(x,2) ~= 1
    x = x';
end

xx = zeros(length(x)+dw-1,1);
x_start = floor(dw/2) + 1;
xx(1:x_start) = x(x_start:-1:1);
xx(x_start+dx:end) = x(dx+x_star-dw:end);
xx(x_start:x_start+dx-1) = x;

for ii = 1:dx
    S(ii,:) = 2*abs(fft(xx(ii:ii+dw-1)))/dw;
end

S = S';

end