% usrpth=userpath;
% usrpth(end)='';
usrpth=fileparts(mfilename('fullpath')); %#ok<GPFST>
cd(usrpth);
addpath(usrpth);
addpath(genpath('C:\\Fiji.app\\scripts'));
warning off MATLAB:dispatcher:nameConflict;
clear all;
close all;
fclose all;
clear java;
clear classes;
if pref('is','Luos'),pref('rm','Luos');end
if pref('is','LuosTrack'),pref('rm','LuosTrack');end
if pref('is','inpdlg'),pref('rm','inpdlg');end
global slash usrpth;
slash=filesep;
pref('set', 'Luos', 'slash', slash); 
usrpth=fileparts(mfilename('fullpath'));
pref('set', 'Luos', 'usrpth', usrpth); 

