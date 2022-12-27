function javaaddpath(varargin)
%JAVAADDPATH Add Java classes to MATLAB
%  JAVAADDPATH DIRJAR adds the specified directory or
%  jar file to the current dynamic Java path. 
%
%  When loading Java classes, MATLAB always searches the static Java path
%  before the dynamic Java path.  The static path is fixed at startup and
%  cannot be changed.  It contains, in the following order:
% 
%     - MATLAB's built-in Java path
%     - the contents of javaclasspath.txt in the startup directory
%     - the contents of javaclasspath.txt in the preferences directory
%       (see 'prefdir')
% 
%  Enter 'javaclasspath' to see the static and current dynamic paths.
%
%  JAVAADDPATH DIRJAR  ... adds directories or jar files
%  to the beginning of the current dynamic Java path.
%  Relative paths are converted to absolute paths.
%
%  JAVAADDPATH ... -END appends the specified directories.
%
%  Use the functional form of JAVAADDPATH, such as 
%  JAVAADDPATH({'dirjar','dirjar',...}), when the directory 
%  specification is stored in a string or cell array of
%  strings.
%
%  Example 1:
%  % Add a directory
%  javaaddpath D:/tools/javastuff 
%  % 'clear java' was used to reload modified Java classes
%
%  Example 2:
%  % Add an internet jar file 
%  javaaddpath http://www.example.com/my.jar
%  javaclasspath % View Java path
%
%  Example 3:
%  % Add the current working directory 
%  javaaddpath(pwd)
%
%  See also JAVACLASSPATH, JAVARMPATH, CLEAR, JAVA. 

% Copyright 2002-2012 The MathWorks, Inc.

n = nargin;

narginchk(1,2);

append = false; % default, pre-append

if n>1
  last = varargin{2};
  
  % append 
  if strcmp(last,'-end')
    append = true;
  
  % pre-append  
  elseif strcmp(last,'-begin')
    append = false;
  end 
end

p = varargin{1};

% save the global vars
a=who('global');
for c=a'
    eval(['global ' c{1} ';']);
    eval(['t.(c{1}) = ' c{1} ';']);
end

% Append or prepend the new path   
if append
  javaclasspath( javaclasspath, p );
else
  javaclasspath( p, javaclasspath );
end

% restore the global vars
for c=a'
    eval(['global ' c{1} ';']);
    eval([c{1} ' = t.(c{1});']);
end

% LocalWords:  DIRJAR dirjar javastuff internet
