import os
import subprocess
from setuptools.command.build_ext import build_ext

class CustomBuild(build_ext):
    def run(self):
        current_dir = os.path.dirname(os.path.abspath(__file__))
        samtools_dir = current_dir
        htslib_dir = os.path.join(samtools_dir, 'htslib-1.9')
        
        # 配置和编译 htslib
        print(f"Configuring htslib in {htslib_dir}")
        subprocess.check_call(['./configure'], cwd=htslib_dir)
        print(f"Compiling htslib in {htslib_dir}")
        subprocess.check_call(['make'], cwd=htslib_dir)
        
        # 配置和编译 samtools
        print(f"Configuring samtools in {samtools_dir}")
        subprocess.check_call(['./configure', f'--with-htslib={htslib_dir}'], cwd=samtools_dir)
        print(f"Compiling samtools in {samtools_dir}")
        subprocess.check_call(['make'], cwd=samtools_dir)
        
        # 安装 samtools
        install_dir = os.path.expanduser('~/biosoft/samtools-1.9')
        print(f"Installing samtools to {install_dir}")
        subprocess.check_call(['make', 'install', f'prefix={install_dir}'], cwd=samtools_dir)
        
        # 运行原始的 build_ext 命令
        build_ext.run(self)