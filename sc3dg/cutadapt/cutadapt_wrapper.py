import os
import sys
import subprocess

def main():
    cutadapt_path = os.path.join(os.path.dirname(__file__), '..', 'bin', 'cutadapt')
    
    if not os.path.exists(cutadapt_path):
        print(f"Error: cutadapt executable not found at {cutadapt_path}", file=sys.stderr)
        sys.exit(1)
    
    command = [cutadapt_path] + sys.argv[1:]
    
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing cutadapt: {e}", file=sys.stderr)
        sys.exit(e.returncode)

if __name__ == "__main__":
    main()