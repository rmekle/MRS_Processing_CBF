%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%disp(sMsg_newLines);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% matExcel_copyValues_btw_Excel_Workbooks_s.m
%
%% Script to copy selected values from different Excel workbooks into another worrkbook
%
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2023;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility 
sFunctionName		= 'matExcel_copyValues_btw_Excel_Workbooks_s';
sMsg_newLines		= sprintf('\n\n');
sMsg_newLine		= sprintf('\n');
disp(sMsg_newLines);


%% Init parameter settings 
% Directory of source files, i.e. of Excel workbooks, from which values are copied from 
% (assuming that all source Excel workbooks are in the same directory)
% Fileneame of one source file
% Directory of filename destination Excel workbook
dirExcelsources			= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_3_2/Diff_Results_Excel_editOFF_water_values/';
fullFileNameExcelsource = '';
dirExceldest			= '/home/mekler/CSB_NeuroRad/mekler/Ralf/CSB_Projects/MRS_Dopamin/MRS_DOPA_DataAnalysis/';
fileNameExceldest		= '3T_MRS_DOPA_Results_2023_Info_Values_mat.xlsx';

% Determine # of source Excel workbooks in directory of source files
source_files			= dir(fullfile(dirExcelsources, '*.xlsx'));
no_source_files			= numel(source_files);

% dir() in MATLAB returns files in an order that can be different from the order from the
% OS, e.g. capital letters are sorted before non-capital letters; hence a sorting
% according to the sorting of the OS is performed
%source_names			= {source_files(:).name}';	% Also possible
source_names			= {source_files.name}';
[~,idx]					= sort(lower(source_names));
%source_names			= source_names(idx)
source_files			= source_files(idx);
warning('%s: Check order of files in list of source files, since ordering in MATLAB can be different from OS!\n', sFunctionName);
disp(sMsg_newLines);

% Specify # of groups, # of value sets (per group), and # of values per vaöue set
% (as a start, # of values per value set is same for all value sets and all groups)
no_groups				= 2;
no_sets					= 2;	% E.g. Mean and STDEV
no_valuesPerSet			= 7;
% no_Sheets				= 3;	% 

% Settings for actually copying all values that were read in and different options for
% source and destination sheets and cell ranges
bCopyValues				= 0;
strSheets_source_groups	= 'Incl1';	% 'All';	% 'Incl1';	% 'Incl2';
strSheets_dest_groups	= 'Incl1';	% 'All';	% 'Incl1';	% 'Incl2';

% Specify sheet(s) and cell ranges in source workbooks, from which values are copied from
% and sheet(s) and cell ranges in destination workbook(s), to which values are copied to
% Source sheets and ranges
% strSheet_source_group_1		= '';
% strSheet_source_group_2		= '',
% strRange_source_group_1		= '';
% strRange_source_group_2		= '';
%acellSheets_source_groups	= cell(no_groups, no_valueSets);
%acellRanges_source_groups	= cell(no_groups, no_valueSets);
switch(strSheets_source_groups)
	case 'All'
		acellSheets_source_groups	= { 'Subjects_Diff_All_L_DOPA', 'Subjects_Diff_All_Placebo'; ...
			'Subjects_Diff_All_L_DOPA', 'Subjects_Diff_All_Placebo'};
		acellRanges_source_groups	= { 'D25:J25', 'D25:J25'; ...
			'D26:J26', 'D26:J26'};
	case 'Incl1'
		acellSheets_source_groups	= { 'Subjects_Diff_Incl1_L_DOPA', 'Subjects_Diff_Incl1_Placebo'; ...
			'Subjects_Diff_Incl1_L_DOPA', 'Subjects_Diff_Incl1_Placebo'};
		acellRanges_source_groups	= { 'D24:J24', 'D24:J24'; ...
			'D25:J25', 'D25:J25'};
	case 'Incl2'
		acellSheets_source_groups	= { 'Subjects_Diff_Incl2_L_DOPA', 'Subjects_Diff_Incl2_Placebo'; ...
			'Subjects_Diff_Incl2_L_DOPA', 'Subjects_Diff_Incl2_Placebo'};
		acellRanges_source_groups	= { 'D23:J23', 'D23:J23'; ...
			'D24:J24', 'D24:J24'};

	otherwise
		error('%s: ERROR: Unknown strSheets_soure_groups =  %s!', sFunctionName, strSheets_source_groups);
end		% End of switch(strSheets_source_groups)

% Destination sheets and ranges
% strSheet_dest_group_1		= '';
% strSheet_dest_group_2		= '';
% strRange_dest_group_1		= '';
% strRange_dest_group_2		= '';
%acellSheets_dest_groups		= cell(no_groups, no_valueSets);
%acellRanges_dest_groups		= cell(no_groups, no_valueSets);
% Hardcoded initialization (for now)
% acellSheets_dest_groups	= { 'MRS_DOPA_3_2_Diff_Conc_All', 'MRS_DOPA_3_2_Diff_Conc_All'; ...
% 							'MRS_DOPA_3_2_Diff_Conc_All', 'MRS_DOPA_3_2_Diff_Conc_All'};
% acellRanges_dest_groups	= { 'B5:H5', 'I5:O5'; ...
% 							'P5:V5', 'W5:AC5'};
% Destination sheets and ranges are selected by a separate switch statement to allow for
% more flexibility; for the L_DOPA study this was not necessary
% For L_DOPA study cell range / starting cellin destination workbook is the same for all 
% options
strRange_dest			= 'B5';
switch(strSheets_dest_groups)
	case 'All'
		strSheet_dest			= 'MRS_DOPA_3_2_Diff_Conc_All';
	case 'Incl1'
		strSheet_dest			= 'MRS_DOPA_3_2_Diff_Conc_Incl1';
	case 'Incl2'
		strSheet_dest			= 'MRS_DOPA_3_2_Diff_Conc_Incl2';

	otherwise
		error('%s: ERROR: Unknown strSheets_dest_groups =  %s!', sFunctionName, strSheets_dest_groups);
end		% End of switch(strSheets_dest_groups)

% Allocate array for all values
% (here: write values for all groups and value sets into same row)
values_all			= zeros(no_source_files, (no_groups*no_sets*no_valuesPerSet));

% For all different source Excel workbooks, read in values for all groups and all value 
% sets into array of all values, reorder array of all values according to order of desired
% output, if required, and then write resulting array (table) of values formatted or 
% using a template into destination workbook
for i=1 : 1 : no_source_files		% no_source_files	% 1
	fprintf('\n\n');
	fprintf('source_files(%d).name = %s\n\n', i, source_files(i).name); 
	fullFileNameExcelsource		= fullfile(dirExcelsources, source_files(i).name);
	%opts_source					= detectImportOptions(fullFileNameExcelsource);
    %preview(fullFileNameExcelsource, opts_source)
	for j=1 : 1 : no_sets
		for k=1 : 1 : no_groups
			fprintf('acellSheets_source_groups{%d, %d} = %s \t acellRanges_source_groups{%d, %d} = %s\n', j, k, acellSheets_source_groups{j, k}, j, k, acellRanges_source_groups{j, k});
			indices					= ( ((j-1)*no_groups*no_valuesPerSet+(k-1)*no_valuesPerSet) + [1:no_valuesPerSet] );
			%opts_source.Sheet		= acellSheets_source_groups{j, k};
			%opts_source.SelectedVariableNames = [4:10];
			%opts_source.DataRange	= 25:25;
			%values_all(i, indices)	= readmatrix(fullFileNameExcelsource, 'Sheet', acellSheets_source_groups{j, k}, 'Range', acellRanges_source_groups{j, k});
			%values_read				= values_all(i, indices)
			%values_read				= readmatrix(fullFileNameExcelsource, 'Sheet', acellSheets_source_groups{j, k}, 'Range', acellRanges_source_groups{j, k})
			values_read				= readmatrix(fullFileNameExcelsource, 'Sheet', acellSheets_source_groups{j, k}, 'Range', acellRanges_source_groups{j, k});
			values_all(i, indices)	= values_read;
			%tmp						= readmatrix(fullFileNameExcelsource, opts_source)
%  			T						= readtable(fullFileNameExcelsource, 'Sheet', acellSheets_source_groups{j, k}, 'Range', acellRanges_source_groups{j, k}, 'ReadVariableNames',false);
%  			disp(T);
% 			M						= readmatrix(fullFileNameExcelsource, 'Sheet', acellSheets_source_groups{j, k}, 'Range', acellRanges_source_groups{j, k});
% 			disp(M);
		end			% End of for k=1 : 1 : no_groups
	end			% End of for j=1 : 1 : no_sets
end		% End of for i=1 : 1 : no_source_files
fprintf('\n\n');
fprintf('values_all = \n');
disp(values_all);

% Write resulting array (table) of values into destination workbook, if selected
fullFileNameExceldest		= fullfile(dirExceldest, fileNameExceldest);
if bCopyValues
	writematrix(values_all, fullFileNameExceldest, 'Sheet', strSheet_dest, 'Range', strRange_dest, 'AutoFitWidth', false);
	%writetable(values_all, fullFileNameExceldest, 'Sheet', strSheet_dest, 'Range', strRange_dest, 'WriteVariableNames', false, 'AutoFitWidth', false);
end		% End of if bCopyValues



%% Code for reading out all different sets of values for each group first and then 
% reordering of values to have values for each set next to each other fro all groups

% % Specify sheet(s) and cell ranges in source workbooks, from which values are copied from
% % and sheet(s) and cell ranges in destination workbook(s), to which values are copied to
% % Source sheets and ranges
% % strSheet_source_group_1		= '';
% % strSheet_source_group_2		= '',
% % strRange_source_group_1		= '';
% % strRange_source_group_2		= '';
% %acellSheets_source_groups	= cell(no_groups, no_valueSets);
% %acellRanges_source_groups	= cell(no_groups, no_valueSets);
% % Hardcoded initialization (for now)
% acellSheets_source_groups	= {'Subjects_Diff_All_L_DOPA', 'Subjects_Diff_All_L_DOPA'; ...
% 								'Subjects_Diff_All_Placebo', 'Subjects_Diff_All_Placebo'};
% acellRanges_source_groups	= {'D25:J25', 'D26:J26'; ...
% 								'D25:J25', 'D26:J26'};
% 
% % Destination sheets and ranges
% % strSheet_dest_group_1		= '';
% % strSheet_dest_group_2		= '';
% % strRange_dest_group_1		= '';
% % strRange_dest_group_2		= '';
% %acellSheets_dest_groups		= cell(no_groups, no_valueSets);
% %acellRanges_dest_groups		= cell(no_groups, no_valueSets);
% % Hardcoded initialization (for now)
% acellSheets_dest_groups	= {'MRS_DOPA_3_2_Diff_Conc_all', 'MRS_DOPA_3_2_Diff_Conc_all'; ...
% 								'MRS_DOPA_3_2_Diff_Conc_all', 'MRS_DOPA_3_2_Diff_Conc_all'};
% acellRanges_dest_groups	= {'B5:H5', 'P5:V5'; ...
% 							'I5:O5', 'W5:AC5'};
% 
% % Allocate array for all values
% % (here: write values for all groups and value sets into same row)
% values_all			= zeros(no_source_files, (no_groups*no_sets*no_valuesPerSet));
% 
% % For all different source Excel workbooks, read in values for all groups and all value 
% % sets into array of all values, reorder array of all values according to order of desired
% % output, and then write resulting array (table) of values formatted or using a template 
% % into destination workbook
% for i=1 : 1 : 2	%no_source_files		% no_source_files	% 1
% 	source_files(i).name
% 	fullFileNameExcelsource		= fullfile(dirExcelsources, source_files(i).name);
% 	for j=1 : 1 : no_groups
% 		for k=1 : 1 : no_sets
% 			acellSheets_source_groups{j, k}
% 			acellRanges_source_groups{j, k}
% 			indices					= ( ((j-1)*no_sets*no_valuesPerSet+(k-1)*no_valuesPerSet) + [1:no_valuesPerSet] );
% 			values_all(i, indices)	= readmatrix(fullFileNameExcelsource, 'Sheet', acellSheets_source_groups{j, k}, 'Range', acellRanges_source_groups{j, k});
% 			values_read				= values_all(i, indices)
% %  			T					= readtable(fullFileNameExcelsource, 'Sheet', acellSheets_source_groups{j, k}, 'Range', acellRanges_source_groups{j, k}, 'ReadVariableNames',false);
% %  			disp(T);
% % 			M					= readmatrix(fullFileNameExcelsource, 'Sheet', acellSheets_source_groups{j, k}, 'Range', acellRanges_source_groups{j, k});
% % 			disp(M);
% 		end			% End of for k=1 : 1 : no_sets
% 	end			% End of for j=1 : 1 : no_groups
% end		% End of for i=1 : 1 : no_source_files
% %disp(values_all);
% 
% % Reorder array of all values according to desired output
% % E.g., mean values for multiple groups should be stored next to each other in each row,
% % followed by values for stdev for all groups;
% % for 2 groups (patients & controls) and 2 sets (mean & stdev), 
% %ind_values			= [1:(no_groups*no_sets*no_valuesPerSet)];
% vec_valuesPerSet	= [1:no_valuesPerSet];
% ind_values_reord	= (repmat(vec_valuesPerSet, (no_groups*no_sets), 1) + ...
% 	[0 (no_sets*no_valuesPerSet) ((no_sets-1)*no_valuesPerSet) ((no_groups-1+no_sets)*no_valuesPerSet)]')';
% ind_values_reord	= ind_values_reord(:)';
% values_all			= values_all(:, ind_values_reord);

