import sys
import subprocess
from pathlib import Path

def cli():
    repo_root = Path(__file__).resolve().parents[1]
    main_py = repo_root / "scripts" / "main.py"
    raise SystemExit(subprocess.call([sys.executable, str(main_py), *sys.argv[1:]]))
