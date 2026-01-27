import os 
from glob import glob 


path_list = glob('./configs/*yaml')

for path in path_list:

	cmd = f'bio-script make -c {path}'

	print(f'Running command: {cmd}')
	os.system(cmd)
	