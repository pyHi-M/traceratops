"""Exercise installer branch selection without network, refs or home writes."""

import os
import subprocess
import tempfile
import unittest
from pathlib import Path


class InstallerTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.bin = self.root / "bin"
        self.bin.mkdir()
        self.parent = self.root / "repositories"
        self.log = self.root / "commands.log"
        self.installer = (
            Path(__file__).resolve().parents[2] / "install_traceratops_uv.bash"
        )
        scripts = {
            "git": """echo "git $*" >> "$TEST_COMMAND_LOG"
case "$1" in
  clone) mkdir -p traceratops ;;
  branch) echo "$TEST_CURRENT_BRANCH" ;;
  status) if [ "$TEST_DIRTY" = yes ]; then echo " M local.py"; fi ;;
esac
""",
            "uv": """echo "uv $*" >> "$TEST_COMMAND_LOG"
if [ "$1" = venv ]; then mkdir -p .venv/bin; touch .venv/bin/activate; fi
""",
            "python": 'echo "python $*" >> "$TEST_COMMAND_LOG"\n',
            "trace_filter": 'echo "trace_filter $*" >> "$TEST_COMMAND_LOG"\n',
        }
        for name, text in scripts.items():
            file = self.bin / name
            file.write_text("#!/usr/bin/env bash\nset -eu\n" + text)
            file.chmod(0o755)

    def run_installer(self, current="dev", dirty="no", branch=None):
        env = dict(os.environ)
        env.update(
            PATH=str(self.bin) + os.pathsep + env["PATH"],
            TRACERATOPS_REPO_DIR=str(self.parent),
            TEST_COMMAND_LOG=str(self.log),
            TEST_CURRENT_BRANCH=current,
            TEST_DIRTY=dirty,
        )
        env.pop("TRACERATOPS_BRANCH", None)
        if branch:
            env["TRACERATOPS_BRANCH"] = branch
        return subprocess.run(
            ["bash", str(self.installer)], env=env, capture_output=True, text=True
        )

    def test_new_clone_selects_dev_and_installs_editable_extra(self):
        result = self.run_installer()
        self.assertEqual(result.returncode, 0, result.stderr)
        log = self.log.read_text()
        self.assertIn(
            "git clone --branch dev https://github.com/pyHi-M/traceratops.git", log
        )
        self.assertIn("uv pip install -e .[dev]", log)

    def test_existing_clean_main_checkout_switches_to_dev(self):
        (self.parent / "traceratops").mkdir(parents=True)
        result = self.run_installer(current="main")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("git switch -- dev", self.log.read_text())
        self.assertNotIn("git pull", self.log.read_text())

    def test_dirty_checkout_refuses_switch_before_replacing_environment(self):
        environment = self.parent / "traceratops/.venv"
        environment.mkdir(parents=True)
        marker = environment / "keep"
        marker.touch()
        result = self.run_installer(current="main", dirty="yes")
        self.assertNotEqual(result.returncode, 0)
        self.assertTrue(marker.exists())
        self.assertNotIn("git switch", self.log.read_text())
        self.assertNotIn("uv venv", self.log.read_text())

    def test_explicit_branch_override(self):
        result = self.run_installer(current="main", branch="main")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("git clone --branch main", self.log.read_text())


if __name__ == "__main__":
    unittest.main()
