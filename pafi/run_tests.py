import sys,os,pytest

def main():
    # Run pytest on the tests directory
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # Construct the path to the tests directory
    tests_dir = os.path.join(current_dir, 'tests')
    exit_code = pytest.main([tests_dir])
    sys.exit(exit_code)

if __name__ == "__main__":
    main()