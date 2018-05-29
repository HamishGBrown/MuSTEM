name = getArgument();
name = replace(name,"(?<!\\\\)\\\\(?!\\\\)", "$0$0");

if (indexOf(name, ".txt") != -1)
{
    name2 = "["+name+"]";
    run("Text Image... ", "open="+name2);  
}
else if (indexOf(name, ".raw") != -1 ^ (indexOf(name, ".bin") !=-1))
{
    bytes=8;
    l=File.length(name);

    fnam = File.getName(name);

	//Check for muSTEM dimensions string in filename
    if(matches(fnam,"^.+_[0-9]+x[0-9]+\..+"))
    {
	//Parse muSTEM dimensions string for dimensions
    x1 =  lastIndexOf(fnam,'_');
    x = lastIndexOf(fnam,'x');
    y = lastIndexOf(fnam,'\.');
    nopiy = parseInt(substring(fnam,x1+1,x));
    nopix = parseInt(substring(fnam,x+1,y));
    }
    else{
	//If no muSTEM dimensions string query user for dimensions
    Dialog.create("Raw input file dimensions");
    Dialog.addNumber("Leading dimension:",sqrt(l/4));
    Dialog.addNumber("Second dimension:",sqrt(l/4));
    Dialog.show();
    nopiy=Dialog.getNumber();
    nopix=Dialog.getNumber();
    }
	
	//Datatype can be intuited from filesize given dimensions
     if (l/(nopiy*nopix) == 1)
    {
    openString="open=["+name+"] image=[8-bit] width="+nopiy+" height="+nopix+" offset=0 number=1 gap=0";
    }
    if (l/(nopiy*nopix) == 2)
    {
    openString="open=["+name+"] image=[16-bit Signed] width="+nopiy+" height="+nopix+" offset=0 number=1 gap=0";
    }
    if (l/(nopiy*nopix) == 8)
    {
    openString="open=["+name+"] image=[64-bit Real] width="+nopiy+" height="+nopix+" offset=0 number=1 gap=0";
    }
    if(l/(nopiy*nopix) == 4)
    {
    openString="open=["+name+"] image=[32-bit Real] width="+nopiy+" height="+nopix+" offset=0 number=1 gap=0";
    }
    run("Raw...",openString);
    
}
else
{
    open(name);
}


