import os
import glob
import sys
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




class fastp_cmd(install):
    def run(self):
        install.run(self)
        self.copy_fastp()

    def copy_fastp(self):
        fastp_source = os.path.join(os.path.dirname(__file__), 'sc3dg', 'utils', 'fastp')
        fastp_dest = os.path.join(self.install_lib, 'sc3dg', 'utils', 'fastp')
        self.copy_file(fastp_source, fastp_dest)
        # 确保 fastp 是可执行的
        os.chmod(fastp_dest, os.stat(fastp_dest).st_mode | stat.S_IEXEC)
 


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

        self.run_command('build_ext')
        self.run_command('build_nanoplexer')
        self.run_command('build_bwa')
        self.run_command('build_bowite2')
        self.run_command('build_minimap2')
        self.run_command('fastp')

        install_samtools()
        install.run(self)







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
    version="0.0.66",
    description="A toolkit for processing single cell Hi-C data",
    author="starker", 
    author_email="caikangwen@126.com", 
    cmdclass={
        'install': CustomInstallCommand,
        'build_ext': bedtools_CustomBuild,
        'build_nanoplexer': NanoplexerCustomBuild,
        'build_bwa': bwa_CustomBuild,
        'build_bowite2': bowtie2_CustomBuild,
    
        'fastp': fastp_cmd,
        'build_minimap2': minimap2_CustomBuild,
        
        

    },
    ext_modules=cythonize([
        Extension(
            'sc3dg.cutadapt.qualtrim',
            sources=['sc3dg/cutadapt/qualtrim.pyx'],
            include_dirs=[np.get_include()],
            define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
        ),

        Extension(
            'sc3dg.cutadapt.info',
            sources=['sc3dg/cutadapt/info.pyx'],
            include_dirs=[np.get_include()],
            define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
        ),
        Extension(
            'sc3dg.cutadapt._align',
            sources=['sc3dg/cutadapt/_align.pyx'],
            include_dirs=[np.get_include()],
            define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
        ),
        Extension(
            'sc3dg.cutadapt._kmer_finder',
            sources=['sc3dg/cutadapt/_kmer_finder.pyx'],
            include_dirs=[np.get_include()],
            define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
        ),
        Extension(
            'sc3dg.model.dyn_util',
            sources=['sc3dg/model/dyn_util.pyx'],
            include_dirs=[np.get_include()],
            define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
        ),
         Extension(
                'sc3dg.pairtools.lib.dedup_cython',
                sources=['sc3dg/pairtools/lib/dedup_cython.pyx'],
                include_dirs=[np.get_include()],
                define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
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