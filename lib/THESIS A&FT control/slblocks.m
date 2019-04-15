function blkStruct = slblocks
%SLBLOCKS Defines the block library for a specific Toolbox or Blockset.

%   Copyright 1986-2002 The MathWorks, Inc. 
% $Revision: 1.13 $

% Name of the subsystem which will show up in the SIMULINK Blocksets
% and Toolboxes subsystem.
% Example:  blkStruct.Name = 'DSP Blockset';
blkStruct.Name = ['Adaptive and ' sprintf('\n') ' Fault Tolerant Control'];

% The function that will be called when the user double-clicks on
% this icon.
% Example:  blkStruct.OpenFcn = 'dsplib';
blkStruct.OpenFcn = 'Adaptive_FaultTolerantControl';%.mdl file

% The argument to be set as the Mask Display for the subsystem.  You
% may comment this line out if no specific mask is desired.
% Example:  blkStruct.MaskDisplay = 'plot([0:2*pi],sin([0:2*pi]));';
blkStruct.MaskDisplay = 'disp(''AFTC'')';

% Define the library list for the Simulink Library browser.
% Return the name of the library model and the name for it
Browser.Library = 'Adaptive_FaultTolerantControl';
Browser.Name    = 'Adaptive and Fault Tolerant Control';
Browser.IsFlat  = 1;

blkStruct.Browser = Browser;
