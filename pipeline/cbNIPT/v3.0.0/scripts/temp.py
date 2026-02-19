import os 


path = '/storage/home/jhkim/scripts/bioscript/generated_scripts'

for p in os.listdir('/storage/home/jhkim/pipelines/cbNIPT/v0.0.3/scripts'):

    x = f'{path}/{p}'

    cmd = f'cp {x} .'
    os.system(cmd)
    print(f'cp {cmd} {cmd.replace()}')