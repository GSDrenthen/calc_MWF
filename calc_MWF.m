%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate Myelin-Water Fraction from 2-D myelin-water imaging data with a 
% slice profile correction using the fourier transform of the excitation pulse
%
% inputs:
%	MWI_data	4-D array of myelin-water imaging data (X,Y,Z,TE)
%	te		1-D vector of echo times, (e.g.: 10:10:320)
%	B1_err		1-D vector of B1 error range, (e.g.: 0.5:0.01:1)
%
% output:
%	MWF		3-D array of myelin-water fraction (X,Y,Z)
%	B1_map		3-D array of B1 error map (X,Y,Z)
%	
% dependencies:
%	SVD_filter.m		(10.1016/j.mri.2006.03.006)
%	calc_sliceprofile.m 	(https://github.com/GSDrenthen/calc_MWF/blob/master/calc_sliceprofile.m)
%  	NonNeg_OMP.m 		(https://github.com/GSDrenthen/Non-Negative-OMP/blob/master/NonNeg_OMP.m)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MWF, B1_map] = calc_MWF(MWI_data, te, B1_err)
	H = fspecial('gaussian',[3 3], 1);
	%H = fspecial3('gaussian',[9 9 9], 1);

	for slice = 1:length(MWI_data(1,1,:,1))
		MWI_data(:,:,slice,:) = SVD_filter(MWI_data(:,:,slice,:));
		for echo = 1:length(MWI_data(1,1,1,:))
			MWI_data(:,:,slice,echo) = imfilter(MWI_data(:,:,slice,echo),H);
		end
	end

	MWI_data = MWI_data + min(MWI_data(:));
	
	T2Times = logspace(log10(te(1)*1.5),log10(2000),1000);
	T2Basis = calc_sliceprofile(te,B1_err,T2Times);

	MWF = zeros(size(MWI_data,1),size(MWI_data,2),size(MWI_data,3));
	B1_map = zeros(size(MWI_data,1),size(MWI_data,2),size(MWI_data,3));

	for zz = 1:size(MWI_data,3)
		for xx = 1:size(MWI_data,1)
			for yy = 1:size(MWI_data,2)
				for B1n = 1:1:length(B1_err)
					[tmp] = lsqnonneg(T2Basis(:,:,B1n),squeeze(MWI_data(xx,yy,zz,:)));
					if size(tmp,1) == 121
						res(B1n,:) = (squeeze(MWI_data(xx,yy,zz,:)) - T2Basis(:,:,B1n)*tmp);
					end
				end  
				[~,B1n] = min(sum(abs(res(:,:)),2));    
				B1_map(xx,yy,zz) = B1n;
			end
		end
	end

	for zz = 1:size(MWI_data,3)
		for xx = 1:size(MWI_data,1)
			for yy = 1:size(MWI_data,2)
				MWF(xx,yy,zz) = NonNeg_OMP(T2Basis(:,:,B1_map(xx,yy,zz)),squeeze(MWI_data(xx,yy,zz,:)),T2Times,20);
			end
		end
	end
end
