import os
import sys
import subprocess

def main():
    # 获取当前脚本所在目录
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 构造minimap2可执行文件的路径
    minimap2_executable = os.path.join(script_dir, 'minimap2')
    
    # 确保minimap2可执行文件是可执行的
    os.chmod(minimap2_executable, 0o755)
    
    # 运行minimap2，传递所有命令行参数
    result = subprocess.run([minimap2_executable] + sys.argv[1:])
    
    # 使用与minimap2相同的退出状态退出
    sys.exit(result.returncode)

if __name__ == '__main__':
    main()