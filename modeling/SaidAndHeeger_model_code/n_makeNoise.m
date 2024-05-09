function [noise] = n_makeNoise(p)

noise = p.noiseamp*randn(1,length(p.tlist));
tkern = normpdf(p.tlist,0,p.noisefilter_t);
tkern = tkern/norm(tkern);
noise = cconv2(noise,tkern);

