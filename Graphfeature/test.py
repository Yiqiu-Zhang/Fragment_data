import subprocess
import os
# export https_proxy="http://10.1.8.5:32680"
local_base = '/mnt/lustre/zhangyiqiu/Fragment_data'
'''
os.system('aws s3 cp s3://Fragment_data/frag/frag_4 /mnt/lustre/zhangyiqiu/Fragment_data/frag/frag_5 --recursive > /dev/null')
'''
with open(os.devnull, 'w') as DEVNULL:
    subprocess.call(['aws', 's3', 'cp', f's3://Fragment_data/frag/frag_4',
                    f'/mnt/lustre/zhangyiqiu/Fragment_data/frag/frag_4','--recursive'],
                    stdout=DEVNULL)#
