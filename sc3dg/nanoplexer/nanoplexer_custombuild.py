import os
import subprocess
from setuptools.command.build_py import build_py

class NanoplexerCustomBuild(build_py):
    def run(self):
        print("NanoplexerCustomBuild: Starting custom build process")
        build_py.run(self)

        nanoplexer_dir = os.path.join(os.path.dirname(__file__))
        print(f"NanoplexerCustomBuild: Nanoplexer directory: {nanoplexer_dir}")
        
        if not os.path.exists(nanoplexer_dir):
            raise FileNotFoundError(f"Nanoplexer directory not found: {nanoplexer_dir}")

        print("NanoplexerCustomBuild: Attempting to compile nanoplexer")
        try:
            subprocess.check_call(['make'], cwd=nanoplexer_dir)
            print("NanoplexerCustomBuild: Compilation successful")
        except subprocess.CalledProcessError as e:
            print(f"Error compiling nanoplexer: {e}")
            raise

        bin_dir = os.path.join(self.build_lib, 'sc3dg', 'bin')
        os.makedirs(bin_dir, exist_ok=True)

        nanoplexer_bin = os.path.join(nanoplexer_dir, 'nanoplexer')
        if os.path.exists(nanoplexer_bin):
            dest = os.path.join(bin_dir, 'nanoplexer')
            self.copy_file(nanoplexer_bin, dest)
            print(f"NanoplexerCustomBuild: Copied nanoplexer to {dest}")
        else:
            print("Warning: nanoplexer binary not found after compilation")

        print("NanoplexerCustomBuild: Custom build process completed")