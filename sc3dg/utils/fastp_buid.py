import os
import shutil
from setuptools import setup, find_packages
from setuptools.command.install import install

class fastp_build(install):
    def run(self):
        install.run(self)
        # 复制 fastp 二进制文件到包的安装目录
        fastp_src = os.path.join(os.path.dirname(__file__), 'sc3dg', 'utils', 'fastp')
        fastp_dst = os.path.join(self.install_lib, 'sc3dg', 'utils', 'fastp')
        os.makedirs(os.path.dirname(fastp_dst), exist_ok=True)
        shutil.copy2(fastp_src, fastp_dst)
        # 确保文件是可执行的
        os.chmod(fastp_dst, 0o755)