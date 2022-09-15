function A = aweight(f)
% computes a weight of input frequencies
% based on https://en.wikipedia.org/wiki/A-weighting


Anum = ((12194^2).*f.^4);

Adiv1 = f.^2 + 20.6^2;
Adiv2 = f.^2 + 107.7^2;
Adiv3 = f.^2 + 737.9^2;
Adiv4 = f.^2 + 12194^2;


A = Anum./(Adiv1.*sqrt(Adiv2.*Adiv3).*Adiv4);
