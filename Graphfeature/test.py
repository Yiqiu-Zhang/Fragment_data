import subprocess
local_base = '/mnt/lustre/zhangyiqiu/Fragment_data'

subprocess.call(['aws', 's3', 'cp', f's3://Fragment_data/feature/feat_{0}',
                 f'{local_base}/feature/feat_{0}/','--recursive'],
                stdout=subprocess.DEVNULL,)