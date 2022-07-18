% function for setting up the common paths across other w3c scripts
function [w3cpath] = fn_w3c_setenvBox()

% set main box directory
w3cpath.BoxMainDir = fn_getBoxDir();

% set path of main w3c script folder
scriptDir = [w3cpath.BoxMainDir, '\YinmingSun\GitHub\w3c-main'];
w3cpath.scriptDir = scriptDir;

% add path to w3c repository subfolders
cd(scriptDir)
addpath('utils','models','data')

end
