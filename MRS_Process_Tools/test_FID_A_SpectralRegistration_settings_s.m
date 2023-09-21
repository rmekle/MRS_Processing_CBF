% FID-A toolkit for MRS
% Test spectral registration settings
%% Geenral parameters
sFunctionname	= 'test-FID_A_SpectralRegistration_settings_s';
dataType		= 'mrs_ref';	% 'mrs_w_ref';	'mrs_ref';	'water_ref';


%% Init parameters used for spectral registration
aaDomain_In				= 'f';			% 'f';		% 't';
tmaxin_In				= 0.2;
iterin_In				= 20;
%alignSS_In				= 2;
aaDomain				= aaDomain_In;
tmaxin					= tmaxin_In;
iterin					= iterin_In;


% Set parameters for drift correction depending on type of data, i.e. whether MRS
% data is spectrum or water signal
% NOTE: Check whether aligning of averages in frequency domain works, if the MR
% spectrum is water signal itself; if not, simply align averages in time domain
%if( strcmp(dataType, 'mrs_w') || strcmp(dataType, 'mrs') )
switch dataType
	case {'mrs', 'mrs_w', 'mrs_w_ref', 'mrs_ref'}
		% MR spectrum is provided together without or with unsuppressed water
		% signal and/or with reference scans
		ppmmin_fix		= 1.6;
		%ppmmaxarray_fix	= [3.5; 4.0; 5.5];
		ppmmaxarray_fix = [2.4,2.85,3.35,4.2,4.4,5.2];
		iamax			= 6;
	case {'water', 'water_ref'}
		% MR spectrum is water signal itself without or with reference scans
		ppmmin_fix		= 4.2;
		ppmmaxarray_fix	= [5.5 5.5 5.2];
		iamax			= 6;

	otherwise
		error('%s: Unknown MRS dataType = %s!', sFunctionName, dataType);
end		% End of switch dataType
noValues_ppmmax		= length(ppmmaxarray_fix);


% now align averages;
%driftCorr=input('Would you like to perform the frequency drift correction?  ','s');
%driftCorr='y';
iter		= 1;
%while (abs(fsPoly(1))>0.001 || abs(phsPoly(1))>0.01) && iter<iterin
while iter<iterin
	fprintf('\n\niter = %d\n', iter);
	%iter			= iter+1
	close all
	%tmax			= 0.25+0.03*randn(1);
	%ppmmin			= 1.6+0.1*randn(1);
	%ppmmaxarray	= [3.5+0.1*randn(1,2),4+0.1*randn(1,3),5.5+0.1*randn(1,1)];
	%ppmmax			= ppmmaxarray(randi(6,1));
	%tmax			= tmaxin+0.03*randn(1)
	tmax			= tmaxin+0.04*randn(1)
	ppmmin			= ppmmin_fix+0.1*randn(1)
	switch noValues_ppmmax
		case 3
			% Generate array of 6 ppmmax values using random number variations
			ppmmaxarray		= [ppmmaxarray_fix(1)+0.1*randn(1,2),ppmmaxarray_fix(2)+0.1*randn(1,3),ppmmaxarray_fix(3)+0.1*randn(1,1)];
		case 6
			% % Generate array of 6 ppmmax values 
			ppmmaxarray		= ppmmaxarray_fix;
			
		otherwise
			error('%s: No option for noValues_ppmmax = %d!', sFunctionName, noValues_ppmmax);
	end		% End of switch noValues_ppmmax
	ppmmax			= ppmmaxarray(randi(iamax,1))
% 	switch aaDomain
% 		case 't'
% 			[out_aa,fs,phs]		= op_alignAverages(out_rm2,tmax,'y');
% 			%[out_aa,fs,phs]		= op_alignAverages(out_rm2);
% 		case 'f'
% 			[out_aa,fs,phs]		= op_alignAverages_fd(out_rm2,ppmmin,ppmmax,tmax,'y');
% 			%[out_aa,fs,phs]		= op_alignAverages_fd(out_rm2,ppmmin,ppmmax,tmax,'n');
% 		otherwise
% 			error('%s: ERROR: avgAlignDomain %s not recognized!', sFunctionName, aaDomain);
% 	end
% 
% 	fsPoly		= polyfit([1:out_aa.sz(out_aa.dims.averages)]',fs,1)
% 	phsPoly		= polyfit([1:out_aa.sz(out_aa.dims.averages)]',phs,1)
% 	%iter
% 	%disp( sprintf('Aligning averages iteration %d', iter) );
% 	fprintf(1, 'Aligning averages iteration %d\n', iter) ;
% 
% 	fscum		= fscum+fs;
% 	phscum		= phscum+phs;
% 
% 	if driftCorr=='y' || driftCorr=='Y'
% 		out_rm2		= out_aa;
% 	end
	iter			= iter+1;
end		% End of while (abs(fsPoly(1))>0.001 || abs(phsPoly(1))>0.01) && iter<iterin