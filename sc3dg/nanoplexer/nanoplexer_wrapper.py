import os
import sys
import subprocess

def main():
    # 获取nanoplexer可执行文件的路径
    nanoplexer_path = os.path.join(os.path.dirname(__file__),  'nanoplexer')
    
    if not os.path.exists(nanoplexer_path):
        print(f"Error: nanoplexer executable not found at {nanoplexer_path}", file=sys.stderr)
        sys.exit(1)
    
    # 构建完整的命令，包括所有参数
    command = [nanoplexer_path] + sys.argv[1:]
    
    try:
        # 执行nanoplexer命令
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing nanoplexer: {e}", file=sys.stderr)
        sys.exit(e.returncode)

if __name__ == "__main__":
    main()