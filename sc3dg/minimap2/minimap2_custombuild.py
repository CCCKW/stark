import os
import subprocess
from setuptools.command.build_ext import build_ext

class CustomBuild(build_ext):
    def run(self):
        # 获取当前目录（setup.py 所在的目录）
        current_dir = os.path.dirname(os.path.abspath(__file__))
        
        # 假设 BWA 源代码就在当前目录下
        minimap2_dir = current_dir
        
        # 编译 BWA
        print(f"Compiling minimap2 in {minimap2_dir}")
        subprocess.check_call(['make', '-C', minimap2_dir])
        
        # 确保 build_lib 目录存在
        if not os.path.exists(self.build_lib):
            os.makedirs(self.build_lib)
        
        # 将编译好的 minimap2 可执行文件复制到 build 目录
        minimap2_executable = os.path.join(minimap2_dir, 'minimap2')
        if os.path.exists(minimap2_executable):
            subprocess.check_call(['cp', minimap2_executable, os.path.join(self.build_lib, 'minimap2')])
        else:
            raise FileNotFoundError(f"Compiled BWA executable not found at {minimap2_executable}")
        
        # 运行原始的 build_ext 命令
        build_ext.run(self)