import re,os,numpy as np

def open_muSTEM_binary(filename):
	'''opens binary with name filename outputted from the muSTEM software'''
	assert os.path.exists(filename), filename+' does not exist'
	#Parse filename for array dimensions
	m = re.search('([0-9]+)x([0-9]+)',filename)
	if m:
		y = int(m.group(2))
		x = int(m.group(1))
	#Get file size and intuit datatype
	size =  os.path.getsize(filename)
	if (size/(y*x) == 4):
		dtype = '>f4'
	elif(size/(y*x) == 8):
		dtype = '>f8'
	print dtype
	#Read data and reshape as required.
	return np.reshape(np.fromfile(filename,dtype = dtype),(y,x))
