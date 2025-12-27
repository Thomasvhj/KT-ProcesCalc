"""
Build executable for Chemical Flow Network GUI

This script creates a standalone .exe file that can be run without Python installed.
Run this script: python build_exe.py
"""

import subprocess
import sys
import os

def main():
    print("=" * 50)
    print("Building Chemical Flow Network GUI Executable")
    print("=" * 50 + "\n")
    
    # Check if PyInstaller is installed
    try:
        import PyInstaller
        print("✓ PyInstaller is installed")
    except ImportError:
        print("Installing PyInstaller...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pyinstaller"])
        print("✓ PyInstaller installed\n")
    
    # Build command
    script_dir = os.path.dirname(os.path.abspath(__file__))
    main_script = os.path.join(script_dir, "chemical_flow_gui.py")
    
    # PyInstaller options:
    # --onefile: Create a single .exe file
    # --windowed: Don't show console window (GUI app)
    # --name: Name of the executable
    # --add-data: Include additional files
    
    cmd = [
        sys.executable, "-m", "PyInstaller",
        "--onefile",
        "--windowed",
        "--name", "ChemicalFlowNetwork",
        "--add-data", f"chemical_flow_network.py;.",
        "--add-data", f"calculations.py;.",
        main_script
    ]
    
    print("Running PyInstaller...")
    print(f"Command: {' '.join(cmd)}\n")
    
    try:
        subprocess.check_call(cmd, cwd=script_dir)
        print("\n" + "=" * 50)
        print("✓ Build successful!")
        print("=" * 50)
        print(f"\nExecutable location:")
        print(f"  {os.path.join(script_dir, 'dist', 'ChemicalFlowNetwork.exe')}")
        print("\nYou can move this .exe file anywhere and run it directly.")
    except subprocess.CalledProcessError as e:
        print(f"\n✗ Build failed: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
