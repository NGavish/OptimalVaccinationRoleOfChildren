%% Sets the appropriate file seperator character (slash or backslash)
%% -----------------------------------------------------------------------------------
%
% Input:
% Non
%
% Output:
%  fileSeperatorChar
% for linux:
% fileSeperatorChar='/'
% for PC:
% fileSeperatorChar='\'
% 
function [commonLibrary] = getCommonLibrary
 if isunix
     commonLibrary='';
 else
     commonLibrary='';
 end
return;


