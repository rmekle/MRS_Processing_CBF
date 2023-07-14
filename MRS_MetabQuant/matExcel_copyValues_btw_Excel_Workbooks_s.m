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
dirExcelsources			= '/home/mekler/CSB_NeuroRad/mekler/Data_II_Analysis/3T_BCAN_MRS_Dopa_Analysis/Z_DOPA_FID-A_SD_3_2/Diff_Results_Excel_editOFF_water/';
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
no_sets		= 2;	% E.g. Mean and STDEV
no_valuesPerSet			= 7;
% no_Sheets				= 3;	% 

% Specify sheet(s) and cell ranges in source workbooks, from which values are copied from
% and sheet(s) and cell ranges in destination workbook(s), to which values are copied to
% Source sheets and ranges
% strSheet_source_group_1		= '';
% strSheet_source_group_2		= '',
% strRange_source_group_1		= '';
% strRange_source_group_2		= '';
%acellSheets_source_groups	= cell(no_groups, no_valueSets);
%acellRanges_source_groups	= cell(no_groups, no_valueSets);
% Hardcoded initialization (for now)
acellSheets_source_groups	= { 'Subjects_Diff_All_L_DOPA', 'Subjects_Diff_All_Placebo'; ...
								'Subjects_Diff_All_L_DOPA', 'Subjects_Diff_All_Placebo'};
acellRanges_source_groups	= { 'D25:J25', 'D25:J25'; ...
								'D26:J26', 'D26:J26'};

% Destination sheets and ranges
% strSheet_dest_group_1		= '';
% strSheet_dest_group_2		= '';
% strRange_dest_group_1		= '';
% strRange_dest_group_2		= '';
%acellSheets_dest_groups		= cell(no_groups, no_valueSets);
%acellRanges_dest_groups		= cell(no_groups, no_valueSets);
% Hardcoded initialization (for now)
acellSheets_dest_groups	= { 'MRS_DOPA_3_2_Diff_Conc_all', 'MRS_DOPA_3_2_Diff_Conc_all'; ...
							'MRS_DOPA_3_2_Diff_Conc_all', 'MRS_DOPA_3_2_Diff_Conc_all'};
acellRanges_dest_groups	= { 'B5:H5', 'I5:O5'; ...
							'P5:V5', 'W5:AC5'};
strSheet_dest			= 'MRS_DOPA_3_2_Diff_Conc_all';
strRange_dest			= 'B5';

% Allocate array for all values
% (here: write values for all groups and value sets into same row)
values_all			= zeros(no_source_files, (no_groups*no_sets*no_valuesPerSet));

% For all different source Excel workbooks, read in values for all groups and all value 
% sets into array of all values, reorder array of all values according to order of desired
% output, if required, and then write resulting array (table) of values formatted or 
% using a template into destination workbook
for i=1 : 1 : 2	%no_source_files		% no_source_files	% 1
	source_files(i).name
	fullFileNameExcelsource		= fullfile(dirExcelsources, source_files(i).name);
	%opts_source					= detectImportOptions(fullFileNameExcelsource);
    %preview(fullFileNameExcelsource, opts_source)
	for j=1 : 1 : no_sets
		for k=1 : 1 : no_groups
			acellSheets_source_groups{j, k}
			acellRanges_source_groups{j, k}
			indices					= ( ((j-1)*no_groups*no_valuesPerSet+(k-1)*no_valuesPerSet) + [1:no_valuesPerSet] );
			%opts_source.Sheet		= acellSheets_source_groups{j, k};
			%opts_source.SelectedVariableNames = [4:10];
			%opts_source.DataRange	= 25:25;
			%values_all(i, indices)	= readmatrix(fullFileNameExcelsource, 'Sheet', acellSheets_source_groups{j, k}, 'Range', acellRanges_source_groups{j, k});
			%values_read				= values_all(i, indices)
			%values_read				= readmatrix(fullFileNameExcelsource, 'Sheet', acellSheets_source_groups{j, k}, 'Range', acellRanges_source_groups{j, k})
			values_read				= readmatrix(fullFileNameExcelsource, 'Sheet', acellSheets_source_groups{j, k}, 'Range', acellRanges_source_groups{j, k})
			values_all(i, indices)	= values_read;
			%tmp						= readmatrix(fullFileNameExcelsource, opts_source)
%  			T						= readtable(fullFileNameExcelsource, 'Sheet', acellSheets_source_groups{j, k}, 'Range', acellRanges_source_groups{j, k}, 'ReadVariableNames',false);
%  			disp(T);
% 			M						= readmatrix(fullFileNameExcelsource, 'Sheet', acellSheets_source_groups{j, k}, 'Range', acellRanges_source_groups{j, k});
% 			disp(M);
		end			% End of for k=1 : 1 : no_groups
	end			% End of for j=1 : 1 : no_sets
end		% End of for i=1 : 1 : no_source_files
disp(values_all);

% Write resulting array (table) of values into destination workbook
fullFileNameExceldest		= fullfile(dirExceldest, fileNameExceldest);
%writematrix(values_all, fullFileNameExceldest, 'Sheet', strSheet_dest, 'Range', strRange_dest, 'AutoFitWidth', false);
%writetable(values_all, fullFileNameExceldest, 'Sheet', strSheet_dest, 'Range', strRange_dest, 'WriteVariableNames', false, 'AutoFitWidth', false);




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

