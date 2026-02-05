from glob import glob
import os


path_list = glob('./configs/*yaml')

for path in path_list:

    cmd = f'bio-script make -c {path} -o ./generated_scripts'
    print(cmd)
    os.system(cmd)