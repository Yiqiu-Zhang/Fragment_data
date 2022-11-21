import subprocess
import os

local_base = '/mnt/lustre/zhangyiqiu/Fragment_data'
'''
os.system('aws s3 cp s3://Fragment_data/frag/frag_4 '
          '/mnt/lustre/zhangyiqiu/Fragment_data/frag/frag_4 --recursive >/dev/null')
'''
subprocess.call(['aws', 's3', 'cp', f's3://Fragment_data/frag/frag_4',
                     f'/mnt/lustre/zhangyiqiu/Fragment_data/frag/frag_4','--recursive', '>/dev/null'],
                    )# stdout=subprocess.DEVNULL
