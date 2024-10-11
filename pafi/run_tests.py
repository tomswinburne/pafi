import sys
import pytest

def main():
    # Run pytest on the tests directory
    exit_code = pytest.main(['tests'])
    sys.exit(exit_code)

if __name__ == "__main__":
    main()