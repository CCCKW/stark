import os
import subprocess
from setuptools.command.build_ext import build_ext

class CustomBuild(build_ext):
    def run(self):
        # 获取当前目录（setup.py 所在的目录）
        current_dir = os.path.dirname(os.path.abspath(__file__))
        
        # 假设 BWA 源代码就在当前目录下
        bedtools_dir = current_dir
        build_output_dir = bedtools_dir + '/bin'
        os.makedirs(build_output_dir, exist_ok=True)

        # 编译 BWA
        print(f"Compiling bedtools in {bedtools_dir}")
        make_command = [
            'make', 
            '-C', bedtools_dir, 
            f'prefix={build_output_dir}'
        ]
        subprocess.check_call(make_command)
        
        # 确保 build_lib 目录存在
        if not os.path.exists(self.build_lib):
            os.makedirs(self.build_lib)
        
        # 将编译好的 bedtools 可执行文件复制到 build 目录
        bedtools_executable = os.path.join(bedtools_dir, 'bin/bedtools')
        if os.path.exists(bedtools_executable):
            subprocess.check_call(['cp', bedtools_executable, os.path.join(self.build_lib, 'bedtools')])
        else:
            raise FileNotFoundError(f"Compiled BWA executable not found at {bedtools_executable}")
        
        # 运行原始的 build_ext 命令
        build_ext.run(self)