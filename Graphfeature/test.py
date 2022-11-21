import subprocess
local_base = '/mnt/lustre/zhangyiqiu/Fragment_data'
for num in range(1,50):
    subprocess.call(['aws', 's3', 'cp', f's3://Fragment_data/frag/frag_4',
                     f'/mnt/lustre/zhangyiqiu/Fragment_data/frag/frag_4','--recursive'],
                    )# stdout=subprocess.DEVNULL
