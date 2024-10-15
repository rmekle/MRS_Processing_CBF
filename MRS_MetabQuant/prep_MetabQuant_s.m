%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% prep_MetabQuant_s.m
%
%% Function to prepare metabolite quantification of magnetic resonance spectroscopy (MRS) data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USAGE
%
%
%



function [status, msg] = prep_MetabQuant_s(outDirString,outDirString_LCM,seqType,options)

%% Parse arguments (required and optional)
arguments
    outDirString		{mustBeText}
	outDirString_LCM	{mustBeText}
    seqType				{mustBeText}
	options.WriteFilenames			(1,1) {islogical}       = 1
	options.CopyFiles				(1,1) {islogical}       = 1
	options.CopyFiles_MRS			(1,1) {islogical}       = 1
	options.CopyFiles_ref_Quant		(1,1) {islogical}       = 0
	options.CopyFiles_ref_ECC		(1,1) {islogical}       = 0
	options.CopyFiles_w				(1,1) {islogical}       = 1
end

% Convert optional inputs to regular variables
bWriteFilenames			= options.WriteFilenames;
bCopyFiles				= options.CopyFiles;
bCopyFiles_MRS			= options.CopyFiles_MRS;
bCopyFiles_ref_Quant	= options.CopyFiles_ref_Quant;
bCopyFiles_ref_ECC		= options.CopyFiles_ref_ECC;
bCopyFiles_w			= options.CopyFiles_w;


%% Set string for name of routine and display blank lines for enhanced output visibility 
% and init output parameters
sFunctionName		= 'prepare_MetabQuant_s';
fprintf('\n\n');
status		= 1;
msg			= '';


%% If output directories for LCM quantification/file copying do not exist, create them
if not(isfolder(outDirString_LCM))
	sMsg = sprintf('%s: Creating output directory %s ...\n', sFunctionName, outDirString_LCM);
    disp(sMsg);
    if ~mkdir(outDirString_LCM)
		error('%s: Could not create (mkdir) output directory %s!\n', sFunctionName, outDirString_LCM);
	end
end
% % Additional output directory for water signals
% if not(isfolder(outDirString_LCM_w))
% 	sMsg = sprintf('%s: Creating output directory %s ...\n', sFunctionName, outDirString_LCM_w);
%     disp(sMsg);
%     if ~mkdir(outDirString_LCM_w)
% 		error('%s: Could not create (mkdir) output directory %s!\n', sFunctionName, outDirString_LCM_w);
% 	end
% end


%% Copy data files depending on sequence type and create list of filenames, if selected
switch seqType
	case 'SPECIAL'
		% Nothing yet
		error('%s: ERROR: No implementation for sequence type %s!', sFunctionName, seqType);
	case 'MEGA-PRESS'
		% Select names of text files that contain lists of files for each signal category
		% MEGA-PRESS editing sequence
		textFileName_diff_MRS	= 'list_filenames_MRS_Diff_Spectra.txt';
		textFileName_OFF		= 'list_filenames_MRS_editOFF_All.txt';
		textFileName_OFF_MRS	= 'list_filenames_MRS_editOFF_Spectra.txt';
		textFileName_OFF_w		= 'list_filenames_MRS_editOFF_Water.txt';

		% Obtain information about the list of different groups of files
		% (assuming that all data files are included in the same directory)
		% (On Linux, file list in Matlab also includes the two directories "." and "..", which
		% means that the actual # of files in the directory is (# of entries in list - 2;
		% however, if dir is used to list specific files, e.g. using a file extension, these two
		% directories are not included in the resulting list)
		%cd(outDirString);
		% MEGA-PRESS difference files (MRS + water)
		structFileListing_diff		= dir([outDirString, '*_diff_*.RAW']);
		noEntriesListing_diff		= length( structFileListing_diff );
		%noDataFiles_diff			= noEntriesListing_diff - 2
		
		% MEGA-PRESS edit_OFF files (MRS + water)
		structFileListing_OFF		= dir([outDirString, '*_editOFF_*.RAW']);
		noEntriesListing_OFF		= length( structFileListing_OFF );
		
		% MEGA-PRESS water difference files (water only)
		structFileListing_diff_w	= dir([outDirString, '*_w_*_diff_*.RAW']);
		noEntriesListing_diff_w		= length( structFileListing_diff_w );
		
		% MEGA-PRESS water edit_OFF files (water only)
		structFileListing_OFF_w		= dir([outDirString, '*_w_*_editOFF_*.RAW']);
		noEntriesListing_OFF_w		= length( structFileListing_OFF_w );
		
		% Determine MEGA-PRESS MRS difference files (MRS only)
		structFileListing_diff_MRS	= structFileListing_diff( ~ismember({structFileListing_diff.name}, {structFileListing_diff_w.name}) );
		noEntriesListing_diff_MRS	= length( structFileListing_diff_MRS );
		
		% Determine MEGA-PRESS MRS edit_OFF files (MRS only)
		structFileListing_OFF_MRS	= structFileListing_OFF( ~ismember({structFileListing_OFF.name}, {structFileListing_OFF_w.name}) );
		noEntriesListing_OFF_MRS	= length( structFileListing_OFF_MRS );

		
		% Copy only MEGA-PRESS MRS difference files and all edit_OFF files with .RAW 
		% extension to destination folder, if selected
		indexStart		= 1;
		indexStep		= 1;
		if bCopyFiles
			% Copy MEGA-PRESS MRS difference files
			for ind=indexStart : indexStep : noEntriesListing_diff_MRS
				% Copy MEGA-PRESS difference file
				[status, msg] = copyfile(fullfile(outDirString, structFileListing_diff_MRS(ind).name), outDirString_LCM);
				if ~status
					disp(msg);
					error('%s: Could not copy file %s to output directory %s!\n', sFunctionName, structFileListing_diff_MRS(ind).name, outDirString_LCM);
				end
			end
			
			% Copy MEGA-PRESS edit_OFF files (MRS + water)
			for ind=indexStart : indexStep : noEntriesListing_OFF
				% Copy MEGA-PRESS edit_OFF file
				[status, msg] = copyfile(fullfile(outDirString, structFileListing_OFF(ind).name), outDirString_LCM);
				if ~status
					disp(msg);
					error('%s: Could not copy file %s to output directory %s!\n', sFunctionName, structFileListing_OFF(ind).name, outDirString_LCM);
				end
			end
			
		end		% End of if bCopyFiles
			
		
		%% Write filenames for specific groups of data files into corresponding textfiles
		% (additional for loops are used to allow for greater flexibility), if selected
		istart		= 1;
		istep		= 1;
		
		if bWriteFilenames
			% Filenames of all MEGA-PRESS edit_OFF files (MRS + water)
			fid_OFF			= fopen(fullfile(outDirString_LCM, textFileName_OFF), 'wt+');
			if fid_OFF == -1
				error('%s: Could not open textfile %s in output directory %s!\n', sFunctionName, textFileName_OFF, outDirString_LCM);
			end
			for i=istart : istep : (noEntriesListing_OFF-1)
				nbytes	= fprintf(fid_OFF, '%s\n', structFileListing_OFF(i).name);
			end
			% Write last filename without new line character at end of line to textfile
			nbytes	= fprintf(fid_OFF, '%s', structFileListing_OFF(noEntriesListing_OFF).name);
			if fclose(fid_OFF) == -1
				error('%s: Could not close textfile %s in output directory %s!\n', sFunctionName, textFileName_OFF, outDirString_LCM);
			end
			
			% Filenames of MEGA-PRESS water edit_OFF files (water only)
			fid_OFF_w		= fopen(fullfile(outDirString_LCM, textFileName_OFF_w), 'wt+');
			if fid_OFF_w == -1
				error('%s: Could not open textfile %s in output directory %s!\n', sFunctionName, textFileName_OFF_w, outDirString_LCM);
			end
			for i=istart : istep : (noEntriesListing_OFF_w-1)
				nbytes	= fprintf(fid_OFF_w, '%s\n', structFileListing_OFF_w(i).name);
			end
			% Write last filename without new line character at end of line to textfile
			nbytes	= fprintf(fid_OFF_w, '%s', structFileListing_OFF_w(noEntriesListing_OFF_w).name);
			if fclose(fid_OFF_w) == -1
				error('%s: Could not close textfile %s in output directory %s!\n', sFunctionName, textFileName_OFF_w, outDirString_LCM);
			end
			
			% Filenames of MEGA-PRESS MRS edit_OFF files (MRS only)
			fid_OFF_MRS		= fopen(fullfile(outDirString_LCM, textFileName_OFF_MRS), 'wt+');
			if fid_OFF_MRS == -1
				error('%s: Could not open textfile %s in output directory %s!\n', sFunctionName, textFileName_OFF_MRS, outDirString_LCM);
			end
			for i=istart : istep : (noEntriesListing_OFF_MRS-1)
				nbytes	= fprintf(fid_OFF_MRS, '%s\n', structFileListing_OFF_MRS(i).name);
			end
			% Write last filename without new line character at end of line to textfile
			nbytes	= fprintf(fid_OFF_MRS, '%s', structFileListing_OFF_MRS(noEntriesListing_OFF_MRS).name);
			if fclose(fid_OFF_MRS) == -1
				error('%s: Could not close textfile %s in output directory %s!\n', sFunctionName, textFileName_OFF_MRS, outDirString_LCM);
			end
			
			% Filenames of MEGA-PRESS MRS difference files (MRS only)
			fid_diff_MRS		= fopen(fullfile(outDirString_LCM, textFileName_diff_MRS), 'wt+');
			if fid_diff_MRS == -1
				error('%s: Could not open textfile %s in output directory %s!\n', sFunctionName, textFileName_diff_MRS, outDirString_LCM);
			end
			for i=istart : istep : (noEntriesListing_diff_MRS-1)
				nbytes	= fprintf(fid_diff_MRS, '%s\n', structFileListing_diff_MRS(i).name);
			end
			% Write last filename without new line character at end of line to textfile
			nbytes	= fprintf(fid_diff_MRS, '%s', structFileListing_diff_MRS(noEntriesListing_diff_MRS).name);
			if fclose(fid_diff_MRS) == -1
				error('%s: Could not close textfile %s in output directory %s!\n', sFunctionName, textFileName_diff_MRS, outDirString_LCM);
			end
			
		end		% End of if bWriteFilenames
	case 'sLASER'
		% Select names of text files that contain lists of files for each signal category
		% sLASER sequence with reference water signals for ECC and ECC
		textFileName_MRS 		= 'list_filenames_MRS_Spectra.txt';
		textFileName_ref_Quant 	= 'list_filenames_MRS_Ref_Quant.txt';
		textFileName_ref_ECC 	= 'list_filenames_MRS_Ref_ECC.txt';
		textFileName_w 			= 'list_filenames_MRS_Water.txt';
		
		% Obtain information about the list of different groups of files
		% (assuming that all data files are included in the same directory)
		
		% sLASER all processed .RAW files
		structFileListing_all		= dir([outDirString, '*_processed_*.RAW']);
		noEntriesListing_all		= length( structFileListing_all );
		%noDataFiles_all			= noEntriesListing_all - 2
		
		% sLASER processed .RAW water reference files for ECC
		structFileListing_ref_ECC	= dir([outDirString, '*_ref_ECC*_processed_*.RAW']);
		noEntriesListing_ref_ECC	= length( structFileListing_ref_ECC );
		
		% sLASER processed .RAW water reference files for quantification (Quant)
		structFileListing_ref_Quant	= dir([outDirString, '*_ref_Quant*_processed_*.RAW']);
		noEntriesListing_ref_Quant	= length( structFileListing_ref_Quant );
				
		% sLASER processed .RAW unsuppressed water files
		structFileListing_w		= dir([outDirString, '*_w_*_processed_*.RAW']);
		noEntriesListing_w		= length( structFileListing_w );
		
 		% sLASER processed .RAW MR spectra files 
		% (find .RAW files that do not belong to any other group)
		structFileListing_Tmp1		= structFileListing_all( ~ismember({structFileListing_all.name}, {structFileListing_ref_Quant.name}) );
		structFileListing_Tmp2		= structFileListing_Tmp1( ~ismember({structFileListing_Tmp1.name}, {structFileListing_ref_ECC.name}) );
 		structFileListing_MRS		= structFileListing_Tmp2( ~ismember({structFileListing_Tmp2.name}, {structFileListing_w.name}) );
 		noEntriesListing_MRS		= length( structFileListing_MRS );

		% Copy (processed) MR Spectra files, water reference signals for quantification,
		% and unsuppresed water signals
		% with . RAW extension to destination folder, if selected
		indexStart		= 1;
		indexStep		= 1;
		if bCopyFiles
			% Copy (processed) MR spectra files, if selected
			if bCopyFiles_MRS
				for ind=indexStart : indexStep : noEntriesListing_MRS
					% Copy MR spectra file
					[status, msg] = copyfile(fullfile(outDirString, structFileListing_MRS(ind).name), outDirString_LCM);
					if ~status
						disp(msg);
						error('%s: Could not copy file %s to output directory %s!\n', sFunctionName, structFileListing_MRS(ind).name, outDirString_LCM);
					end
				end
			end		% End of if bCopyFiles_MRS
			
			% Copy (processed) water reference files quantification (Quant), if selected
			if bCopyFiles_ref_Quant
				for ind=indexStart : indexStep : noEntriesListing_ref_Quant
					% Copy water reference file quantification (Quant)
					[status, msg] = copyfile(fullfile(outDirString, structFileListing_ref_Quant(ind).name), outDirString_LCM);
					if ~status
						disp(msg);
						error('%s: Could not copy file %s to output directory %s!\n', sFunctionName, structFileListing_ref_Quant(ind).name, outDirString_LCM);
					end
				end
			end		% End of if bCopyFiles_ref_Quant
			
			% Copy (processed) water reference files for ECC, if selected
			if bCopyFiles_ref_ECC
				for ind=indexStart : indexStep : noEntriesListing_ref_ECC
					% Copy water reference file for ECC
					[status, msg] = copyfile(fullfile(outDirString, structFileListing_ref_ECC(ind).name), outDirString_LCM);
					if ~status
						disp(msg);
						error('%s: Could not copy file %s to output directory %s!\n', sFunctionName, structFileListing_ref_ECC(ind).name, outDirString_LCM);
					end
				end
			end		% End of if bCopyFiles_ref_ECC
			
			% Copy (processed) unsuppressed water files, if selected
			if bCopyFiles_w
				for ind=indexStart : indexStep : noEntriesListing_w
					% Copy unsuppressed water file
					[status, msg] = copyfile(fullfile(outDirString, structFileListing_w(ind).name), outDirString_LCM);
					if ~status
						disp(msg);
						error('%s: Could not copy file %s to output directory %s!\n', sFunctionName, structFileListing_w(ind).name, outDirString_LCM);
					end
				end
			end		% End of if bCopyFiles_w
			
		end		% End of if bCopyFiles
		
		
		%% Write filenames for specific groups of data files into corresponding textfiles
		% (additional for loops are used to allow for greater flexibility), if selected
		istart		= 1;
		istep		= 1;
		if bWriteFilenames
			% Filenames of (processed) MR spectra files
			if bCopyFiles_MRS
				fid_MRS			= fopen(fullfile(outDirString_LCM, textFileName_MRS), 'wt+');
				if fid_MRS == -1
					error('%s: Could not open textfile %s in output directory %s!\n', sFunctionName, textFileName_MRS, outDirString_LCM);
				end
				for i=istart : istep : (noEntriesListing_MRS-1)
					nbytes	= fprintf(fid_MRS, '%s\n', structFileListing_MRS(i).name);
				end
				% Write last filename without new line character at end of line to textfile
				nbytes	= fprintf(fid_MRS, '%s', structFileListing_MRS(noEntriesListing_MRS).name);
				if fclose(fid_MRS) == -1
					error('%s: Could not close textfile %s in output directory %s!\n', sFunctionName, textFileName_MRS, outDirString_LCM);
				end
			end		% End of if bCopyFiles_MRS
			
			% Filenames of (processed) water reference files for quantification (Quant)
			if bCopyFiles_ref_Quant
				fid_ref_Quant	= fopen(fullfile(outDirString_LCM, textFileName_ref_Quant), 'wt+');
				if fid_ref_Quant == -1
					error('%s: Could not open textfile %s in output directory %s!\n', sFunctionName, textFileName_ref_Quant, outDirString_LCM);
				end
				for i=istart : istep : (noEntriesListing_ref_Quant-1)
					nbytes	= fprintf(fid_ref_Quant, '%s\n', structFileListing_ref_Quant(i).name);
				end
				% Write last filename without new line character at end of line to textfile
				nbytes	= fprintf(fid_ref_Quant, '%s', structFileListing_ref_Quant(noEntriesListing_ref_Quant).name);
				if fclose(fid_ref_Quant) == -1
					error('%s: Could not close textfile %s in output directory %s!\n', sFunctionName, textFileName_ref_Quant, outDirString_LCM);
				end
			end		% End of if bCopyFiles_ref_Quant
			
			% Filenames of (processed) water reference files for ECC
			if bCopyFiles_ref_ECC
				fid_ref_ECC		= fopen(fullfile(outDirString_LCM, textFileName_ref_ECC), 'wt+');
				if fid_ref_ECC == -1
					error('%s: Could not open textfile %s in output directory %s!\n', sFunctionName, textFileName_ref_ECC, outDirString_LCM);
				end
				for i=istart : istep : (noEntriesListing_ref_ECC-1)
					nbytes	= fprintf(fid_ref_ECC, '%s\n', structFileListing_ref_ECC(i).name);
				end
				% Write last filename without new line character at end of line to textfile
				nbytes	= fprintf(fid_ref_ECC, '%s', structFileListing_ref_ECC(noEntriesListing_ref_ECC).name);
				if fclose(fid_ref_ECC) == -1
					error('%s: Could not close textfile %s in output directory %s!\n', sFunctionName, textFileName_ref_ECC, outDirString_LCM);
				end
			end		% End of if bCopyFiles_ref_ECC
			
			% Filenames of (processed) unsuppressed water files
			if bCopyFiles_w
				fid_w			= fopen(fullfile(outDirString_LCM, textFileName_w), 'wt+');
				if fid_w == -1
					error('%s: Could not open textfile %s in output directory %s!\n', sFunctionName, textFileName_w, outDirString_LCM);
				end
				for i=istart : istep : (noEntriesListing_w-1)
					nbytes	= fprintf(fid_w, '%s\n', structFileListing_w(i).name);
				end
				% Write last filename without new line character at end of line to textfile
				nbytes	= fprintf(fid_w, '%s', structFileListing_w(noEntriesListing_w).name);
				if fclose(fid_w) == -1
					error('%s: Could not close textfile %s in output directory %s!\n', sFunctionName, textFileName_w, outDirString_LCM);
				end
			end		% End of if bCopyFiles_w
			
		end		% End of if bWriteFilenames
		
	otherwise
		error('%s: ERROR: Unknown sequence type %s!', sFunctionName, seqType);
end		% End of switch seqType


