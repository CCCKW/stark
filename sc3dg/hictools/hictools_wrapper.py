#!/usr/bin/env python3
import os
import sys
import subprocess

def main():
    # 获取 hictools 可执行文件的路径
    current_dir = os.path.dirname(os.path.abspath(__file__))
    hictools_path = os.path.join(current_dir, 'hictools')
    
    if not os.path.exists(hictools_path):
        print(f"Error: hictools executable not found at {hictools_path}")
        sys.exit(1)
    
    # 执行 hictools 并传递所有参数
    cmd = [hictools_path] + sys.argv[1:]
    subprocess.run(cmd)

if __name__ == '__main__':
    main()