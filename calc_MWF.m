%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Myelin-Water Fraction from mutli-echo GRASE data
%
% dependencies:
%  SVD_filter.m
%  RegNNLS.m
%  calc_sliceprofile.m
% 

function [MWF, B1_map] = calc_MWF(MWI_data, te, B1_err)
	addpath('dependencies')

	% MWI_data -> X,Y,Z,TE

	H = fspecial('gaussian',[9 9], 1);
	%H = fspecial3('gaussian',[9 9 9], 1);

	for slice = 1:length(MWI_data(1,1,:,1))
		MWI_data(:,:,slice,:) = SVD_filter(MWI_data(:,:,slice,:));
		for echo = 1:length(MWI_data(1,1,1,:))
			MWI_data(:,:,slice,echo) = imfilter(MWI_data(:,:,slice,echo),H);
		end
	end

	MWI_data = MWI_data + min(MWI_data(:));
	
	T2Times = logspace(log10(te(1)*1.5),log10(2000),120);
	%T2Times = logspace(log10(te(1)*1.5),log10(2000),1000);
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
				x = RegNNLS(T2Basis(:,:,B1_map(xx,yy,zz)),squeeze(MWI_data(xx,yy,zz,:)),1.02,1.025);
				MWF(xx,yy,zz) = sum(x(T2Times<40)) / sum(x);
				MWF(xx,yy,zz) = NonNeg_OMP(T2Basis(:,:,B1_map(xx,yy,zz)),squeeze(MWI_data(xx,yy,zz,:)),T2Times,10);
			end
		end
	end
end