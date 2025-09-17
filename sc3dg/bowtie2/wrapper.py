import os
import sys
import subprocess

def main():
    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Construct the path to the Perl script
    perl_script = os.path.join(script_dir, 'bowtie2')
    
    # Make sure the Perl script is executable
    os.chmod(perl_script, 0o755)
    
    # Run the Perl script with all command-line arguments
    result = subprocess.run(['perl', perl_script] + sys.argv[1:])
    
    # Exit with the same status as the Perl script
    sys.exit(result.returncode)

if __name__ == '__main__':
    main()