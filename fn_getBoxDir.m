% Gets the path for the main Box directory.
% May need to be configured based on setup on each computer

function [BoxMainDir] = fn_getBoxDir()

if ismac
    % Mac OS: get current username
    username = char(java.lang.System.getProperty('user.name'));
    
    % BoxMainDir = ['/Users/', username, '/Box']; %Mac OS settings 1
    BoxMainDir = ['/Users/', username, '/Library/CloudStorage/Box-Box/']; %Mac OS settings 2
    
elseif ispc
    
    username=getenv('USERNAME');
    BoxMainDir = ['C:\Users\', username, '\Box\'];
    
end


end