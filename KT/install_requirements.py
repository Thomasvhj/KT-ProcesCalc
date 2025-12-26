"""
Install required packages for Chemical Flow Network GUI

Run this script to ensure all dependencies are installed:
    python install_requirements.py
"""

import subprocess
import sys

def install_package(package):
    """Install a package using pip"""
    print(f"Installing {package}...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])
    print(f"✓ {package} installed successfully\n")

def check_package(package, import_name=None):
    """Check if a package is installed"""
    if import_name is None:
        import_name = package
    try:
        __import__(import_name)
        return True
    except ImportError:
        return False

def main():
    print("=" * 50)
    print("Chemical Flow Network GUI - Dependency Installer")
    print("=" * 50 + "\n")
    
    # Required packages: (pip_name, import_name)
    required = [
        ("numpy", "numpy"),
    ]
    
    # Check Python version
    print(f"Python version: {sys.version}")
    if sys.version_info < (3, 7):
        print("⚠ Warning: Python 3.7+ is recommended for dataclasses support")
    print()
    
    # Check and install packages
    all_installed = True
    for pip_name, import_name in required:
        if check_package(pip_name, import_name):
            print(f"✓ {pip_name} is already installed")
        else:
            all_installed = False
            try:
                install_package(pip_name)
            except subprocess.CalledProcessError:
                print(f"✗ Failed to install {pip_name}")
                print(f"  Try manually: pip install {pip_name}")
    
    # Check tkinter (comes with Python but not always)
    print("\nChecking tkinter...")
    try:
        import tkinter
        print("✓ tkinter is available")
    except ImportError:
        print("✗ tkinter is NOT available")
        print("  tkinter usually comes with Python.")
        print("  On Linux, try: sudo apt-get install python3-tk")
        print("  On macOS, reinstall Python from python.org")
        all_installed = False
    
    print("\n" + "=" * 50)
    if all_installed:
        print("All dependencies are installed!")
        print("Run the GUI with: python chemical_flow_gui.py")
    else:
        print("Some dependencies may need manual installation.")
    print("=" * 50)

if __name__ == "__main__":
    main()
