import os
import sys
import subprocess

def main():
    fastp_path = os.path.join(os.path.dirname(__file__), 'fastp')
    if not os.path.exists(fastp_path):
        print(f"Error: fastp executable not found at {fastp_path}", file=sys.stderr)
        sys.exit(1)
    
    try:
        subprocess.run([fastp_path] + sys.argv[1:], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running fastp: {e}", file=sys.stderr)
        sys.exit(e.returncode)

if __name__ == '__main__':
    main()