function [ mE_interp, mG_interp, tau_interp ] = interpolateTLMparams( request_beta, request_lambda, k_request )
%[ mE_interp, mG_interp, tau_interp ] = interpolateTLMparams( request_beta, request_lambda, k )
%   This interpolates TLM transfer function parameters based on beta
%   (dissipation number), lambda (taper ratio), and number of parameters k.
%   Refer to ven der Buhs and Wiens, 2017 for details.
%
%   References
%   J ven der Buhs and T Wiens. Modelling Dynamic Response of Hydraulic Fluid Within Tapered Transmission Lines. Proceedings of the 15th Scandinavian International Conference on Fluid Power, 2017

load TaperedTLMLookupTable.mat

if k_request ~= 6
    error('Currenly only 6 weighting factors are available, k must be 6 for this version of the model')
end

mE_interp=nan(1,k_request);
mG_interp=nan(1,k_request);

for i=1:k_request
    mE_interp(i) = interp2(beta_lookup,lambda_lookup,mE_lookup(:,:,i),request_beta,request_lambda,'linear');
    mG_interp(i) = interp2(beta_lookup,lambda_lookup,mG_lookup(:,:,i),request_beta,request_lambda,'linear');
end

tau_interp = interp2(beta_lookup,lambda_lookup,tau_lookup,request_beta,request_lambda,'linear');
end

