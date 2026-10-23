import importlib.util
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PATHS_FILE = ROOT / "paths.py"


def load_paths():
    specification = importlib.util.spec_from_file_location("ultrapipe_paths", PATHS_FILE)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def test_repository_path_defaults_to_checkout(monkeypatch):
    monkeypatch.delenv("ULTRAPIPE_PATH", raising=False)
    paths = load_paths()

    assert os.environ["ULTRAPIPE_PATH"] == str(ROOT)
    assert paths.get_repository_path() == ROOT
    assert paths.get_data_path() == ROOT / "data"
    assert paths.get_visuals_path() == ROOT / "visuals"


def test_repository_path_honors_environment_override(tmp_path, monkeypatch):
    monkeypatch.setenv("ULTRAPIPE_PATH", str(tmp_path))

    assert load_paths().get_repository_path() == tmp_path.resolve()