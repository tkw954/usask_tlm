function [ t11, t12, t21, t22 ] = TaperedTLMTransferMatrix( omega, E_1, E_2, F_1, F_2, G_1, G_2, Zc_1, Zc_2, T_1, T_2 )
%[ t11, t12, t21, t22 ] = TaperedTLMTransferMatrix( omega, E_1, E_2, F_1, F_2, G_1, G_2, Zc_1, Zc_2, T_1, T_2 )
% Calculates the transmission matrix from the TLM transfer functions,
% characteristic impedance, and wave propagation times.
%
% Inputs: Transfer functions, characteristic impedance, and wave
%          propagation times.
%
% Outputs: t11,t21,t12,t22 = TLM transmission matrix terms
%
% NOTE: E, F, G, T, and Z_c are separated by _1 and _2 to allow investigation into
% separate transfer functions as discussed in ven der Buhs and Wiens [2017].
%
% References:
% J ven der Buhs and T Wiens. Modelling Dynamic Response of Hydraulic Fluid Within Tapered Transmission Lines. 
% Proceedings of the 15th Scandinavian International Conference on Fluid Power, 2017
%
% N Johnston, M Pan, and S Kudzma. An enhanced transmission line method for modelling laminar flow of liquid 
% in pipelines. Journal of Systems and Control Engineering, 228(4):193–206, 2014.

t11 = ((E_1+Zc_1).*G_1.^(-1).*exp(1j*omega*T_1)+F_1.*G_2.*exp(-1j*omega*T_2))./(E_1+Zc_1+F_1); %eq (20)

t21 = ((-G_2.*exp(-1j.*omega.*T_2)+G_1.^(-1).*exp(1j*omega*T_1)))./(E_1+Zc_1+F_1); %eq (22)

t12 = ((E_1+Zc_1).*(E_2+Zc_2).*G_1.^(-1).*exp(1j*omega*T_1)-F_1.*F_2.*G_2.*exp(-1j*omega*T_2))./((E_1+Zc_1+F_1)); %eq (21)

t22 = -((E_2+Zc_2).*G_1.^(-1).*exp(1j*omega*T_1)+F_2.*G_2.*exp(-1j*omega*T_2))./(E_1+Zc_1+F_1); %eq (23)

end

