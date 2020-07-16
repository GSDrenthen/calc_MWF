%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate the T2 basis set (T2 decay profiles for varying T2 times)
%
% inputs:
%	te		1-D vector of echo times, (e.g.: 10:10:320)
%	B1_err		1-D vector of B1 error range, (e.g.: 0.5:0.01:1)
%	T2Times		1-D vector of T2 times, (e.g.: logspace(log10(te(1)*1.5),log10(2000),1000))
%
% outputs
%	T2Basis		3-D array of T2 decay profiles (TE, T2 times, B1 errors)
%
% dependencies:
%	cp_cpmg_epg_domain_fplus_fminus.m (10.1002/jmri.24619)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function T2Basis = calc_sliceprofile(te,B1_err,T2Times)
	Sinc_Gauss = [66; -138; -359; -595; -845; -1107; -1378; -1658; -1943; -2230;
	-2517; -2800; -3076; -3341; -3592; -3824; -4034; -4217; -4370; -4488;
	-4567; -4602; -4591; -4530; -4414; -4240; -4007; -3710; -3348; -2918;
	-2421; -1853; -1216; -510; 266; 1110; 2019; 2991; 4024; 5114;
	6256; 7446; 8678; 9947; 11247; 12571; 13911; 15262; 16614; 17960;
	19293; 20604; 21885; 23129; 24326; 25471; 26554; 27570; 28511; 29371;
	30144; 30826; 31410; 31894; 32274; 32547; 32712; 32767; 32712; 32547;
	32274; 31894; 31410; 30826; 30144; 29371; 28511; 27570; 26554; 25471;
	24326; 23129; 21885; 20604; 19293; 17960; 16614; 15262; 13911; 12571;
	11247; 9947; 8678; 7446; 6256; 5114; 4024; 2991; 2019; 1110;
	266];
	num_sample = 101;
	t = linspace(0,4.224e-3,num_sample);
	max_val = 4.01;
	pulse_dur = 4.224e-3;
	Sinc_Gauss = (max_val .* Sinc_Gauss) ./ (max(Sinc_Gauss));

	Fs = num_sample / pulse_dur;
	T = 1/Fs;
	L = num_sample;
	t = (0:L-1)*T;
	y = Sinc_Gauss;
	n = 2^nextpow2(L);
	n = 1000;
	Y = fft(y,n)/L;
	f = Fs*(0:(n/2))/n;
	P = abs(Y/n);

	sp = [P(ceil(360.5/f(2)):-1:1); P(1:ceil(360.5/f(2)))];
	sp = sp / max(sp);

	T2Basis_sp = zeros(32,length(sp));
	T2Basis = zeros(32,length(T2Times),length(B1_err));

	for T2 = 1:1:length(T2Times)
		for ii = 1:length(B1_err)
			Mex = sin(deg2rad(90 * sp * B1_err(ii)));
			for nn = 1:1:length(sp)
				T2Basis_sp(:,nn) = Mex(nn) * (cp_cpmg_epg_domain_fplus_fminus(length(te),[180; 150.*ones(31,1)].*sp(nn).*B1_err(ii) ,te(1),1000,T2Times(T2)));
			end
			T2Basis(:,T2,ii) = sum(T2Basis_sp,2);
		end
	end
end
