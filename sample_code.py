def setting_arg():
	parser = argparse.ArgumentParser(description='FlowGrid Parameters')
	parser.add_argument('--f',type=str,required=True,dest="file",
		help='input the file name')
	parser.add_argument('--n',type=int,required=True,dest="bin_n",
		help='input the number of bins')
	parser.add_argument('--eps',type=float,required=True,dest="eps",
		help='input the maximun distance for connected bins')
	parser.add_argument('--t',type=int,required=False,dest="MinDenC",
		help='input the minimuns collective density for core bins')
	parser.add_argument('--o',type=str,required=False,dest="output",
		help='input the output location')
	parser.add_argument('--l',type=str,required=False,dest="label",
		help='input the label location')
	pars=parser.parse_args()
	file=pars.file
	bin_n=pars.bin_n
	MinDenC=pars.MinDenC
	eps=pars.eps
	output=pars.output
	label=pars.label
	return file,bin_n,MinDenC,eps,output,label

def check_file_valid(file):
	if file[-4:]!=".csv":
		print("FlowGrid code only accept csv file.")
		os._exit()
	if not  os.path.isfile(file):
		print(file+" does not exist, please input the correct location of file")
		os._exit()
	data=np.genfromtxt(file, delimiter=',',skip_header=1)
	if data.shape[1]==1:
		print("The format of "+file+" may not be wrong.")
		os._exit()
	print("The number of cells is: "+ str(data.shape[0]))
	print("The number of dimensions is: "+ str(data.shape[1]))
	return data

if __name__ == "__main__":
	from FlowGrid import *
	import numpy as np
	from time import time
	import argparse
	import os
	file,bin_n,MinDenC,eps,output,label_file=setting_arg()
	data=check_file_valid(file)
	t1=time()
	if MinDenC:
		fg=FlowGrid(data,bin_n=bin_n,eps=eps, MinDenC=MinDenC)
	else:
		fg=FlowGrid(data,bin_n=bin_n,eps=eps)
	label=fg.clustering()
	print("runing time: "+ str(round(time()-t1,3)))
	if output:
		np.savetxt(output, label, delimiter=',',fmt="%d")
	else:
		np.savetxt(file[:-4]+"_FlowGrid_label.csv", label, delimiter=',',fmt="%d")
	if label_file:
		from sklearn.metrics import adjusted_rand_score as ARI
		true_label=np.genfromtxt(label_file, delimiter=',')
		print("ARI:"+ str(round(ARI(true_label,label),4)))