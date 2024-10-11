import subprocess
import sys

def check_lammps():
    try:
        # Try to import lammps
        import lammps
        print("LAMMPS is installed and accessible.")
        return True
    except ImportError:
        print("LAMMPS is not installed or not accessible.")
        return False

def suggest_lammps_installation():
    print("\nTo install LAMMPS, please use conda:")
    print("conda install -c conda-forge lammps")
    print("\nAfter installing LAMMPS, you can install pafi using:")
    print("pip install pafi")

def main():
    if not check_lammps():
        suggest_lammps_installation()
        sys.exit(1)
    
    print("All dependencies are satisfied. You can proceed with using pafi.")

if __name__ == "__main__":
    main()