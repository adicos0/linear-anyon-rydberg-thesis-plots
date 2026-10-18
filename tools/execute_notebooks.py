"""Execute default notebook cells in separate kernels with bounded timeouts.

Run from the repository root with the environment containing the requirements.
Kernel argv is bound to the running Python, avoiding accidental global kernels.
"""
from pathlib import Path
import os
import sys
from tempfile import TemporaryDirectory
import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager

ROOT = Path(__file__).resolve().parents[1]


class EnvironmentKernelManager(KernelManager):
    def format_kernel_cmd(self, extra_arguments=None):
        return [sys.executable, "-m", "ipykernel_launcher", "-f", self.connection_file,
                *(extra_arguments or [])]


def main():
    paths = sorted((ROOT / "notebooks").glob("*.ipynb"))
    if len(sys.argv) > 1:
        selected = set(sys.argv[1:])
        unknown = selected - {path.name for path in paths}
        if unknown:
            raise SystemExit(f"Unknown notebook filenames: {sorted(unknown)}")
        paths = [path for path in paths if path.name in selected]
    # Isolate runtime caches from both the repository and the user's settings.
    with TemporaryDirectory(prefix="thesis-notebooks-") as cache:
        for key, directory in (("IPYTHONDIR", "ipython"),
                               ("MPLCONFIGDIR", "matplotlib"),
                               ("XDG_CACHE_HOME", "cache")):
            os.environ.setdefault(key, str(Path(cache) / directory))
        execute_all(paths)


def execute_all(paths):
    for path in paths:
        notebook = nbformat.read(path, as_version=4)
        nbformat.validate(notebook)
        client = NotebookClient(notebook, timeout=180, kernel_name="python3",
                                kernel_manager_class=EnvironmentKernelManager,
                                resources={"metadata": {"path": str(ROOT)}},
                                allow_errors=False)
        client.execute()
        nbformat.validate(notebook)
        # Do not expose local runtime paths in the committed record.
        for cell in notebook.cells:
            cell.get("metadata", {}).pop("execution", None)
        nbformat.write(notebook, path)
        count = sum(cell.cell_type == "code" for cell in notebook.cells)
        print(f"PASS {path.name}: {count} executed code cells", flush=True)


if __name__ == "__main__":
    main()
