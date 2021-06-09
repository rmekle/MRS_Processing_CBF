% op_average_subSpecs_s.m
% Jamie Near, McGill University 2014, modified by Ralf Mekle, Charite Berlin 2021.
% 
% USAGE:
% out=op_average_subSpecs_s(in);
% 
% DESCRIPTION:
% Combine the subspecs in a scan by adding the subspecs together and then 
% dividing by the number of subspecs.
% 
% INPUTS:
% in	= input data in matlab structure format.
%
% OUTPUTS:
% out   = Output following averaging.  

function out = op_average_subSpecs_s(in)

%if in.flags.averaged || in.subspecs<2
if in.subspecs<2
    %DO NOTHING
    disp('WARNING: No subSpecs found (in.subspecs < 2). Returning input without modification!');
    return;
end

if in.dims.subSpecs==0
    %DO NOTHING
    disp('WARNING: No subSpecs found (in.dims.subSpecs = 0). Returning input without modification!');
    out=in;
    return;
else
    
    % add the spectra along the subSpecs dimension;
    fids=sum(in.fids,in.dims.subSpecs);
    fids=squeeze(fids);
    fids=fids/in.sz(in.dims.subSpecs); %divide by number of subSpecs;
    
    %re-calculate Specs using fftveraged
    specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);
    
    %change the dims variables.
    if in.dims.t>in.dims.subSpecs
        dims.t=in.dims.t-1;
    else
        dims.t=in.dims.t;
    end
    if in.dims.coils>in.dims.subSpecs
        dims.coils=in.dims.coils-1;
    else
        dims.coils=in.dims.coils;
    end
    dims.subSpecs=0;
    if in.dims.averages>in.dims.subSpecs
        dims.averages=in.dims.averages-1;
    else
        dims.averages=in.dims.averages;
    end
    if in.dims.extras>in.dims.subSpecs
        dims.extras=in.dims.extras-1;
    else
        dims.extras=in.dims.extras;
    end
    
    
    %re-calculate the sz variable
    sz=size(fids);
    
    
    %FILLING IN DATA STRUCTURE
    out=in;
    out.fids=fids;
    out.specs=specs;
    out.sz=sz;
    out.dims=dims;
    out.subspecs=1;
	out.averages=in.averages/in.subspecs;	% Assuming that in.averages is multiple of 
											% in.subspecs
											
    %FILLING IN THE FLAGS
    out.flags=in.flags;
    out.flags.writtentostruct=1;
    out.flags.isISIS=0;						% No more ISIS data available after averaging

end



