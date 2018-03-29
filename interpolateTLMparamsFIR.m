function [ mE_interp, mG_interp, tau_interp, K_diff, K_mean ] = interpolateTLMparamsFIR( request_beta, request_lambda, k_request )
%[ mE_interp, mG_interp, tau_interp ] = interpolateTLMparams( request_beta, request_lambda, k )
%   This interpolates TLM transfer function parameters based on beta
%   (dissipation number), lambda (taper ratio), and number of parameters k.
%   Refer to ven der Buhs and Wiens, 2017 for details.
%
%   References
%   J ven der Buhs and T Wiens. Modelling Dynamic Response of Hydraulic Fluid Within Tapered Transmission Lines. Proceedings of the 15th Scandinavian International Conference on Fluid Power, 2017

if nargin<3
    k_request=6;
end


tmp=load('FIRdata.mat');
FIRdata=tmp.FIRdata;

if k_request ~= 6
    error('Currenly only 6 weighting factors are available, k must be 6 for this version of the model')
end

mE_interp=nan(1,k_request);
mG_interp=nan(1,k_request);

for i=1:k_request
    mE_interp(i) = interp2(FIRdata.lambda,FIRdata.beta,FIRdata.mE(:,:,i),request_lambda,request_beta,'linear');
    mG_interp(i) = interp2(FIRdata.lambda,FIRdata.beta,FIRdata.mG(:,:,i),request_lambda,request_beta,'linear');
end

tau_interp = interp2(FIRdata.lambda,FIRdata.beta,FIRdata.tau,request_lambda,request_beta,'linear');
K_diff = interp2(FIRdata.lambda,FIRdata.beta,FIRdata.K_diff,request_lambda,request_beta,'linear');
K_mean = interp2(FIRdata.lambda,FIRdata.beta,FIRdata.K_mean,request_lambda,request_beta,'linear');

end

