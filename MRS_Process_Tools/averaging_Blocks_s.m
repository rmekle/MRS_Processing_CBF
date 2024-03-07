%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% averaging_Blocks_s.m
% Ralf Mekle, Charite Universitätsmedizin Berlin, Germany, 2024;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% USAGE:
% out = averaging_Blocks_s(in, blockSize);
% 
% DESCRIPTION:
% Average blocks of averages in a scan by adding the averages together and then 
% dividing by the number of averages in one block
% 
% INPUTS:
% in			= Input data in matlab structure format
% blockSize		= Size of blocks (i.e. # of averages in one block) used for averaging
%
% OUTPUTS:
% out			= Output in matlab structure format after averaging of blocks of data 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function out = averaging_Blocks_s(in, blockSize)

%% Clear all variables from workspace and close all figures
% clear all;
% close all;


%% Set string for name of routine and display blank lines for enhanced output visibility
sFunctionName		= 'averaging_Blocks_s';


%% Check whether averaging is possible
if in.flags.averaged || in.averages<2 || in.dims.averages==0
	% Averaging of blocks not possible
    % DO NOTHING
    %disp('WARNING: No averages found. Returning input without modification!');
	warning('%s: No averages found. Returning input without modification!\n\n', sFunctionName);
    out		= in;
    return;
end


% Check whether total # of averages in input data is a multiple of the block size and
% issue warning, if it is not
noAveragesIn	= in.sz(in.dims.averages);
if mod(noAveragesIn, blockSize) ~= 0
	warning('%s: No of averages = %d found is not a multiple of block size = %d. Last block will have fewer averages!\n\n', ...
		sFunctionName, noAveragesIn, blockSize);
end


%% Average block size of averages and concatenate averaged blocks to obtain one output 
% data structure

% If block size is equal to 1, just perform regular averaging
if blockSize == 1
	out		= op_averaging(in);
else
	% For total # of blocks, extract # of averages equal to block size, average each 
	% block, and concatenate blocks together
	out				= struct([]);
	out_block		= struct([]);
	out_block_avg	= struct([]);
	noBlocks		= ceil(noAveragesIn/blockSize);
	for i=1 : 1 : noBlocks
		% Determine indices of each subsequent block of averages and take into account
		% that the maximum index of a block is equal to total # of averages, if total # of
		% averages is not a multiple of block size
		maxIndex		= min(i*blockSize, noAveragesIn);
		indicesBlock	= [((i-1)*blockSize + 1):maxIndex];
		
		% Extract indexed averages for each block from original data, average the averages
		% of each block, and concatenate these averaged blocks into output data structure
		out_block		= op_takeaverages(in, indicesBlock);
		out_block_avg	= op_averaging(out_block);
		out				= op_concatAverages(out, out_block_avg);
	end			% End of for i=1 : 1 : noBlocks

end		% End of if blockSize == 1


