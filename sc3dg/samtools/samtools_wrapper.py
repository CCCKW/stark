#!/usr/bin/env python
import os
import sys
import subprocess

def main():
    # 获取安装的samtools路径
    # /cluster/home/Kangwen/anaconda3/envs/stark/lib/python3.9/site-packages/sc3dg/samtools/bin/samtools
    # /cluster/home/Kangwen/Hic/AAA_tutorial/sc3dg/sc3dg/samtools
    samtools_path = os.path.join(os.path.dirname(__file__), 'bin', 'samtools')
    
    if not os.path.exists(samtools_path):
        print(f"Error: Samtools executable not found at {samtools_path}", file=sys.stderr)
        sys.exit(1)
    
    # 构建完整的命令，包括所有参数
    command = [samtools_path] + sys.argv[1:]
    
    try:
        # 执行samtools命令
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing samtools: {e}", file=sys.stderr)
        sys.exit(e.returncode)

if __name__ == "__main__":
    main()