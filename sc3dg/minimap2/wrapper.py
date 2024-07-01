import os
import subprocess
import sys

def main():
    # 获取包内 minimap2 可执行文件的路径
    minimap2_path = os.path.join(os.path.dirname(__file__), 'minimap2', 'minimap2')
    
    # 确保 minimap2 可执行文件存在
    if not os.path.exists(minimap2_path):
        print("Error: minimap2 executable not found.", file=sys.stderr)
        sys.exit(1)
    
    # 运行 minimap2 命令,传递所有命令行参数
    try:
        subprocess.run([minimap2_path] + sys.argv[1:], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running minimap2: {e}", file=sys.stderr)
        sys.exit(e.returncode)

if __name__ == "__main__":
    main()