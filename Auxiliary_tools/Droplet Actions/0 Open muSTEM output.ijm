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
    else if(matches(fnam,"[0-9]+x[0-9]+_.+"))
    {
    x1 =  indexOf(fnam,'_');
    x = indexOf(fnam,'x');
    nopiy = parseInt(substring(fnam,0,x));
    nopix = parseInt(substring(fnam,x+1,x1));
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
     bytes =  (l-8)/(nopiy*nopix);   
    offset=bytes;
    openString = "width="+nopiy+" height="+nopix+" offset="+offset+" number=1 gap=0";
 
     if (bytes  == 1)
    {
    openString="open=["+name+"] image=[8-bit] "+openString ;
    }
    if (bytes == 2)
    {
    openString="open=["+name+"] image=[16-bit Signed] "+openString ;
    }
    if (bytes == 8)
    {
    openString="open=["+name+"] image=[64-bit Real] "+openString ;
    }
    if(bytes  == 4)
    {
    openString="open=["+name+"] image=[32-bit Real] "+openString ;
    }
    
}
else
{
    open(name);
}


