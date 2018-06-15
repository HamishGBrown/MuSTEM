function [ array_out ] = read_muSTEM_output( title,fmt )
% The function binary_in(title,fmt) reads in a binary file outputed from the MuSTEM
% software. The variable keyword should be a string containing the filepath. The variable
% fmt is optional and specifies the endianness of the file. Use 'b' for Big-endian (default) and 'l'
% for Little-endian
    %Get filename as seperate from directory
	[~,name,~] = fileparts(title);
     
	%Assert the existence of the file
    assert(exist(title,'file')~=1,strcat('Could not find file ',title));
	
	%Pass filename for array dimensions
    dimensions = regexp(name,'([0-9]+)x([0-9]+)','tokens');
    dimensions = dimensions{1};
    Npx = str2num(dimensions{1});
    Npy = str2num(dimensions{2});
   
    %With knowledge of the array dimensions, the file datatype can be intuited   
    s = dir(title);
    type = s.bytes/Npx/Npy;
    if(type==8)
        format='double';
    else
        format='single';
    end
   
    %Get endianness
    if (nargin>1)
        fmt_= fmt;
    else
        fmt_ = 'b';
    end
   
    %Read array from binary file
    fileID = fopen(title);
    array_out = fread(fileID,[Npx Npy],format,0,fmt_);
    fclose(fileID);
end