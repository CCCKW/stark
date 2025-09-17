import os
import glob
import sys
import stat
import shutil
import numpy as np
import subprocess
from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools.command.egg_info import egg_info
from setuptools.command.build_py import build_py
from setuptools.command.build_ext import build_ext
from sc3dg.bwa.custom_build import CustomBuild as bwa_CustomBuild    
from sc3dg.bowtie2.bowtie2_custombulid import CustomBuild as bowtie2_CustomBuild
from sc3dg.minimap2.minimap2_custombuild import CustomBuild as minimap2_CustomBuild
from sc3dg.nanoplexer.nanoplexer_custombuild import NanoplexerCustomBuild
from sc3dg.bedtools.bedtools_build import CustomBuild as bedtools_CustomBuild 
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
# chmod u+x -R sc3dg
cmd = 'chmod u+x -R sc3dg'
os.system(cmd)
os.makedirs('sc3dg/bedtools/bin', exist_ok=True)





def install_samtools():
    samtools_dir = os.path.join(os.path.dirname(__file__), 'sc3dg', 'samtools')
    print(f"Samtools directory: {samtools_dir}")
    if not os.path.exists(samtools_dir):
        raise FileNotFoundError(f"Samtools directory not found: {samtools_dir}")
    
    # 获取 Conda 环境的路径
    conda_env_path = os.environ.get('CONDA_PREFIX', os.path.dirname(os.path.dirname(__file__)))
    install_dir = os.path.join(conda_env_path, 'lib', 'python3.9', 'site-packages', 'sc3dg', 'samtools')
    os.makedirs(install_dir, exist_ok=True)
    
    home_dir = os.path.expanduser('~')
    bzip2_include = os.path.join(home_dir, 'opt', 'bzip2', 'include')
    bzip2_lib = os.path.join(home_dir, 'opt', 'bzip2', 'lib')
    
    configure_command = [
        './configure',
        f'CPPFLAGS=-I{bzip2_include}',
        f'LDFLAGS=-L{bzip2_lib} -Wl,-R{bzip2_lib}',
        '--without-curses',
        '--disable-bz2',
        '--disable-lzma',
        f'--prefix={install_dir}'
    ]
    
    print("Running configure command...")
    subprocess.check_call(configure_command, cwd=samtools_dir)
    print("Running make command...")
    subprocess.check_call(['make'], cwd=samtools_dir)
    print("Running make install command...")
    subprocess.check_call(['make', 'install'], cwd=samtools_dir)
    print(f"Samtools installed in {install_dir}")

    


class CustomInstallCommand(install):
    def run(self):
        # 编译并安装 BWA
        self.build_and_install_bwa()
        
   
        
        # 编译并安装 minimap2
        self.build_and_install_minimap2()
        
        # 编译并安装 Nanoplexer
        self.build_and_install_nanoplexer()
        self.build_and_install_bedtools()


        # 编译并安装 hictools
        self.build_and_install_hictools()
        # 安装 samtools
        install_samtools()
        
        
        # 运行原始的安装命令
        install.run(self)
        
        # 复制 fastp
        self.copy_fastp()


    def build_and_install_hictools(self):
        """编译和安装 hictools"""
        print("Starting hictools build process")
        hictools_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sc3dg', 'hictools')
        
        if not os.path.exists(hictools_dir):
            raise FileNotFoundError(f"hictools directory not found: {hictools_dir}")
        
        # 检查源文件是否存在
        hictools_cpp = os.path.join(hictools_dir, 'hictools.cpp')
        if not os.path.exists(hictools_cpp):
            raise FileNotFoundError(f"hictools.cpp not found in {hictools_dir}")
        
        print(f"Compiling hictools in {hictools_dir}")
        try:
            # 编译 hictools
            compile_command = ['g++', '-g', 'hictools.cpp', '-o', 'hictools']
            subprocess.check_call(compile_command, cwd=hictools_dir)
            print("hictools compilation successful")
        except subprocess.CalledProcessError as e:
            print(f"Error compiling hictools: {e}")
            raise
        
        # 复制编译好的可执行文件
        hictools_executable = os.path.join(hictools_dir, 'hictools')
        if os.path.exists(hictools_executable):
            hictools_dest = os.path.join(self.install_lib, 'sc3dg', 'hictools', 'hictools')
            os.makedirs(os.path.dirname(hictools_dest), exist_ok=True)
            shutil.copy2(hictools_executable, hictools_dest)
            print(f"Copied hictools from {hictools_executable} to {hictools_dest}")
            
            # 设置执行权限
            current_mode = os.stat(hictools_dest).st_mode
            os.chmod(hictools_dest, current_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
            print(f"Set execute permissions for {hictools_dest}")
        else:
            raise FileNotFoundError(f"Compiled hictools executable not found at {hictools_executable}")
        
    
    
    
    
    
    
    
    
    
    def build_and_install_bedtools(self):
        print("Starting bedtools build process")
        bedtools_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sc3dg', 'bedtools')
        build_output_dir = os.path.join(bedtools_dir, 'bin')
        os.makedirs(build_output_dir, exist_ok=True)

        print(f"Compiling bedtools in {bedtools_dir}")
        make_command = [
            'make', 
            '-C', bedtools_dir, 
            f'prefix={build_output_dir}'
        ]
        try:
            subprocess.check_call(make_command)
            print("bedtools compilation successful")
        except subprocess.CalledProcessError as e:
            print(f"Error compiling bedtools: {e}")
            raise

        bedtools_executable = os.path.join(build_output_dir, 'bedtools')
        if os.path.exists(bedtools_executable):
            bedtools_dest = os.path.join(self.install_lib, 'sc3dg', 'bedtools', 'bedtools')
            os.makedirs(os.path.dirname(bedtools_dest), exist_ok=True)
            shutil.copy2(bedtools_executable, bedtools_dest)
            print(f"Copied bedtools from {bedtools_executable} to {bedtools_dest}")
            
            current_mode = os.stat(bedtools_dest).st_mode
            os.chmod(bedtools_dest, current_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
            print(f"Set execute permissions for {bedtools_dest}")
        else:
            raise FileNotFoundError(f"Compiled bedtools executable not found at {bedtools_executable}")
            
    def build_and_install_bwa(self):
        # BWA 编译和安装的代码（与之前相同）
        current_dir = os.path.dirname(os.path.abspath(__file__))
        bwa_dir = os.path.join(current_dir, 'sc3dg', 'bwa')
        
        print(f"Compiling BWA in {bwa_dir}")
        subprocess.check_call(['make', '-C', bwa_dir])
        
        bwa_executable = os.path.join(bwa_dir, 'bwa')
        if os.path.exists(bwa_executable):
            bwa_dest = os.path.join(self.install_lib, 'sc3dg', 'bwa', 'bwa')
            os.makedirs(os.path.dirname(bwa_dest), exist_ok=True)
            shutil.copy2(bwa_executable, bwa_dest)
            print(f"Copied BWA from {bwa_executable} to {bwa_dest}")
            
            current_mode = os.stat(bwa_dest).st_mode
            os.chmod(bwa_dest, current_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
            print(f"Set execute permissions for {bwa_dest}")
        else:
            raise FileNotFoundError(f"Compiled BWA executable not found at {bwa_executable}")


    def build_and_install_bowtie2(self):
        current_dir = os.path.dirname(os.path.abspath(__file__))
        bowtie2_dir = os.path.join(current_dir, 'sc3dg', 'bowtie2')
        
        print(f"Compiling Bowtie2 in {bowtie2_dir}")
        subprocess.check_call(['make', '-C', bowtie2_dir])
        
        bowtie2_executables = ['bowtie2', 'bowtie2-build', 'bowtie2-inspect']
        for executable in bowtie2_executables:
            source = os.path.join(bowtie2_dir, executable)
            if os.path.exists(source):
                dest = os.path.join(self.install_lib, 'sc3dg', 'bowtie2', executable)
                os.makedirs(os.path.dirname(dest), exist_ok=True)
                shutil.copy2(source, dest)
                print(f"Copied Bowtie2 executable from {source} to {dest}")
                
                current_mode = os.stat(dest).st_mode
                os.chmod(dest, current_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
                print(f"Set execute permissions for {dest}")
            else:
                raise FileNotFoundError(f"Compiled Bowtie2 executable not found at {source}")


    def build_and_install_minimap2(self):
        current_dir = os.path.dirname(os.path.abspath(__file__))
        minimap2_dir = os.path.join(current_dir, 'sc3dg', 'minimap2')
        
        print(f"Compiling minimap2 in {minimap2_dir}")
        subprocess.check_call(['make', '-C', minimap2_dir])
        
        minimap2_executable = os.path.join(minimap2_dir, 'minimap2')
        if os.path.exists(minimap2_executable):
            minimap2_dest = os.path.join(self.install_lib, 'sc3dg', 'minimap2', 'minimap2')
            os.makedirs(os.path.dirname(minimap2_dest), exist_ok=True)
            shutil.copy2(minimap2_executable, minimap2_dest)
            print(f"Copied minimap2 from {minimap2_executable} to {minimap2_dest}")
            
            current_mode = os.stat(minimap2_dest).st_mode
            os.chmod(minimap2_dest, current_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
            print(f"Set execute permissions for {minimap2_dest}")
        else:
            raise FileNotFoundError(f"Compiled minimap2 executable not found at {minimap2_executable}")


    def build_and_install_nanoplexer(self):
        print("Starting Nanoplexer build process")
        nanoplexer_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sc3dg', 'nanoplexer')
        print(f"Nanoplexer directory: {nanoplexer_dir}")
        
        if not os.path.exists(nanoplexer_dir):
            raise FileNotFoundError(f"Nanoplexer directory not found: {nanoplexer_dir}")

        print("Attempting to compile Nanoplexer")
        try:
            subprocess.check_call(['make'], cwd=nanoplexer_dir)
            print("Nanoplexer compilation successful")
        except subprocess.CalledProcessError as e:
            print(f"Error compiling Nanoplexer: {e}")
            raise

        nanoplexer_bin = os.path.join(nanoplexer_dir, 'nanoplexer')
        if os.path.exists(nanoplexer_bin):
            nanoplexer_dest = os.path.join(self.install_lib, 'sc3dg', 'nanoplexer', 'nanoplexer')
            os.makedirs(os.path.dirname(nanoplexer_dest), exist_ok=True)
            shutil.copy2(nanoplexer_bin, nanoplexer_dest)
            print(f"Copied Nanoplexer from {nanoplexer_bin} to {nanoplexer_dest}")
            
            current_mode = os.stat(nanoplexer_dest).st_mode
            os.chmod(nanoplexer_dest, current_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
            print(f"Set execute permissions for {nanoplexer_dest}")
        else:
            raise FileNotFoundError(f"Compiled Nanoplexer executable not found at {nanoplexer_bin}")

    def copy_fastp(self):
        # fastp 复制的代码（与之前相同）
        fastp_source = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'sc3dg', 'utils', 'fastp')
        
        if self.root:
            fastp_dest = os.path.join(self.root, self.install_lib[1:], 'sc3dg', 'utils', 'fastp')
        else:
            fastp_dest = os.path.join(self.install_lib, 'sc3dg', 'utils', 'fastp')
        
        os.makedirs(os.path.dirname(fastp_dest), exist_ok=True)
        
        try:
            shutil.copy2(fastp_source, fastp_dest)
            print(f"Copied fastp from {fastp_source} to {fastp_dest}")
            
            current_mode = os.stat(fastp_dest).st_mode
            os.chmod(fastp_dest, current_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
            print(f"Set execute permissions for {fastp_dest}")
        except IOError as e:
            print(f"Error copying fastp: {e}")
            raise
        except OSError as e:
            print(f"Error setting permissions for fastp: {e}")
            raise

    def get_outputs(self):
        outputs = super().get_outputs()
        # 添加 BWA, Bowtie2, minimap2 和 fastp 到输出列表
        bwa_dest = os.path.join(self.install_lib, 'sc3dg', 'bwa', 'bwa')
        bowtie2_dest = os.path.join(self.install_lib, 'sc3dg', 'bowtie2', 'bowtie2')
        minimap2_dest = os.path.join(self.install_lib, 'sc3dg', 'minimap2', 'minimap2')
        fastp_dest = os.path.join(self.install_lib, 'sc3dg', 'utils', 'fastp')
        for dest in [bwa_dest, bowtie2_dest, minimap2_dest, fastp_dest]:
            if os.path.exists(dest):
                outputs.append(dest)
        return outputs







# 获取NumPy包含目录
numpy_includes = np.get_include()

import importlib.util
spec = importlib.util.spec_from_file_location(
    "sc3dg.samtools.htslib-1.19.1",
    "sc3dg/samtools/htslib-1.19.1/__init__.py"
)
htslib = importlib.util.module_from_spec(spec)
spec.loader.exec_module(htslib)


with open('requirements.txt') as f:
    install_requires = f.read().splitlines()


def get_ext_modules():
    extensions = []
    
    
    parse_pysam_path = os.path.join("sc3dg", "pairtools", "lib", "parse_pysam.pyx")
    if os.path.exists(parse_pysam_path):
        import pysam
        extensions.append(
            Extension(
                "sc3dg.pairtools.lib.parse_pysam",
                [parse_pysam_path],
                extra_link_args=pysam.get_libraries(),
                include_dirs=pysam.get_include() + [numpy_includes],
                define_macros=pysam.get_defines() + [('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
                language="c",
                 extra_compile_args=['-std=c99'],
            )
        )
    

    src_files = glob.glob(
        os.path.join(os.path.dirname(__file__), "sc3dg", "pairtools", "lib", "*.pyx")
    )

    for src_file in src_files:
        name = "sc3dg.pairtools.lib." + os.path.splitext(os.path.basename(src_file))[0]
        if "pysam" not in name and "regions" not in name and name != "sc3dg.pairtools.lib.parse_pysam":
            extensions.append(Extension(name, [src_file], 
                                        include_dirs=[numpy_includes],
                                         extra_compile_args=['-std=c99'],
                                        define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')]))
        elif "regions" in name:
            extensions.append(
                Extension(
                    name,
                    [src_file],
                    language="c++",
                    include_dirs=[numpy_includes],
                    define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
                )
            )

    return extensions






setup(
    name="sc3dg",
    version="1.0.0",
    description="A toolkit for processing single cell Hi-C data",
    author="starker", 
    author_email="caikangwen@126.com", 
    cmdclass={
        'install': CustomInstallCommand,
      
        
        

    },
    ext_modules=cythonize([
        Extension(
            'sc3dg.cutadapt.qualtrim',
            sources=['sc3dg/cutadapt/qualtrim.pyx'],
            include_dirs=[np.get_include()],
            define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
            extra_compile_args=['-std=c99'],

        ),

        Extension(
            'sc3dg.cutadapt.info',
            sources=['sc3dg/cutadapt/info.pyx'],
            include_dirs=[np.get_include()],
            define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
            extra_compile_args=['-std=c99'],

        ),
        Extension(
            'sc3dg.cutadapt._align',
            sources=['sc3dg/cutadapt/_align.pyx'],
            include_dirs=[np.get_include()],
            define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
            extra_compile_args=['-std=c99'],

        ),
        Extension(
            'sc3dg.cutadapt._kmer_finder',
            sources=['sc3dg/cutadapt/_kmer_finder.pyx'],
            include_dirs=[np.get_include()],
            define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
            extra_compile_args=['-std=c99'],

        ),
        Extension(
            'sc3dg.model.dyn_util',
            sources=['sc3dg/model/dyn_util.pyx'],
            include_dirs=[np.get_include()],
            define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
            extra_compile_args=['-std=c99'],

        ),
         Extension(
                'sc3dg.pairtools.lib.dedup_cython',
                sources=['sc3dg/pairtools/lib/dedup_cython.pyx'],
                include_dirs=[np.get_include()],
                define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
                extra_compile_args=['-std=c99'],

            ),
             Extension(
                'sc3dg.pairtools.lib.regions',
                sources=['sc3dg/pairtools/lib/regions.pyx'],
                include_dirs=[np.get_include()],
                extra_compile_args=['-std=c++11'],
                language="c++",
                
            ),
         Extension(
        "sc3dg.cooltools.cooltools.lib._numutils",
        ["sc3dg/cooltools/cooltools/lib/_numutils.pyx"],
        include_dirs=[np.get_include()],
        define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
        extra_compile_args=['-std=c99'],

    ),
    ] + get_ext_modules()) ,
    packages=find_packages(exclude=['sc3dg.samtools.htslib-1.19.1']) +
                            ['sc3dg.cutadapt', 'sc3dg.nanoplexer']+['sc3dg.cutadapt']+
                         ['sc3dg', 'sc3dg.utils', 'sc3dg.analysis',  
                            'sc3dg.commands', 
                            'sc3dg.model',
                            'sc3dg.pairtools',
                            'sc3dg.bowtie2', 
                            'sc3dg.minimap2',
                            'sc3dg.bwa', 
                            'sc3dg.cooltools',
                             'sc3dg.cooltools.cooltools.lib',
                             'sc3dg.cooltools.cooltools.api',
                             'sc3dg.cooltools.cooltools.sandbox',
                             'sc3dg.samtools',
                             'sc3dg.bedtools',
                             'sc3dg.bedtools.bin'
                        
        
                             ],
    package_data={
        'sc3dg.bedtools': ['bin/*'],
        'sc3dg.bedtools': ['**/*'], 
        'sc3dg.cutadapt': ['*.py', '*.pyi', '*.so'],
        'sc3dg.nanoplexer': ['*'],
       'sc3dg.cutadapt': ['*'],
        'sc3dg': [
            'bedtool/bin/*',
            'utils/*',
            'utils/fastp',
            'utils/BarcodeIdentification_v1.2.0.jar',
            'model/*.pyx', 'pairtools/lib/*.pyx',
            'bowtie2/*', 'minimap2/*', 'bwa/*',
            'cooltools/*',
            'cooltools/cooltools/api/*.pyx',
            'cooltools/cooltools/cli/*',
            'cooltools/cooltools/lib/*',
            'cooltools/cooltools/*',
            'sc3dg/samtools/*',
            'bismark/*'],
        'sc3dg.cooltools.cooltools.lib': ['*.pyx', '*.pxd'],
        'sc3dg.samtools': ['htslib-1.19.1/**/*'],

        'sc3dg.bismark': ['bismark/**/*'],
        'sc3dg.bowtie2': ['bowtie2/**/*'],
        'sc3dg.minimap2': ['minimap2/**/*'],
        'sc3dg.bwa': ['bwa/**/*'],
        
    },
    install_requires=install_requires,
    setup_requires=['numpy', 'cython'],
    entry_points={
        'console_scripts': [
            'stark = sc3dg.cli:cli',
            'pairtools = sc3dg.pairtools.cli:cli',
            'bowtie2 = sc3dg.bowtie2.wrapper:main',
            'bwa = sc3dg.bwa.bwa_wrapper:main',
            'minimap2 = sc3dg.minimap2.minimap2_wrapper:main',
            'nanoplexer = sc3dg.nanoplexer.nanoplexer_wrapper:main',
            'bismark = sc3dg.bismark.bismark_wrapper:main',
            'cooltools = sc3dg.cooltools.cli:cli',
            'samtools = sc3dg.samtools.samtools_wrapper:main',
            'fastp = sc3dg.utils.fastp_wrapper:main',
            'cutadapt = sc3dg.cutadapt.cli:main_cli',
            'bedtools = sc3dg.bedtools.bedtools_wrapper:main',

        ]
    },
)