dataType		= 'mrs_w';
filename		= filename_svs_slaser_dkd_mrs_w;
filenamew		= filename_w_svs_slaser_dkd_mrs_w;
dirString		= [dirStringBase dataType '/'];
outDirString	= dirString;


% Optional, if water unsuppressed signal exists
strOVS_w		= 'woutOVS';		% 'wOVS';


% Optional, if OVS setting for acquisition of MR spectrum was different
strOVS			= 'woutOVS';		% 'wOVS';


% Optional, if sequence type for acquisition of MR spectrum was different
seqType			= 'PRESS';			% 'sLASER';


% Optional, if any of the following parameters should be changed
aaDomain		= 't';			% 'f';
iterin			= 20;			% Any number probably up to 100
nSD				= 4.0;			% 3.2;		% 2.6;
plotSwitch		= 1;			% 0;
reportSwitch	= 0;			% 1;
strMinUserIn	= 'n';			% 'y';
tmaxin			= 0.2;			% 






