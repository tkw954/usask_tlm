function [ni, k2, mgi, mei, deltai] = gettlmcoeffs(b, nterms)
% Function to get parameters for TLM model for a specific dissipation
% number b
%
% (c) D N Johnston, University of Bath 
% 25 July 2013
%
if nterms < 3 || nterms > 7
    error 'Invalid number of terms'
end

load tlm_weights
%Perform interpolation on sqrt(beta) for slightly better accuracy.
sbeta = sqrt(tlmdata{nterms}.beta);
sb = sqrt(b);

me = tlmdata{nterms}.me;
mg = tlmdata{nterms}.mg;
delta = tlmdata{nterms}.delta;

ni = 0.3*[3.^[0:nterms-1]]/(1+b*3);

mei = interp1(sbeta, me, sb, 'linear', 'extrap');
mgi = interp1(sbeta, mg, sb, 'linear', 'extrap');
deltai = interp1(sbeta, delta, sb, 'linear', 'extrap');
k2 = 1 - 8*b / sum(mei./ni);
