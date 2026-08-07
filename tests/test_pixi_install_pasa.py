import os
import shutil
import stat
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PIX_INSTALL_PASA = REPO_ROOT / "install_scripts" / "pixi_install_pasa.sh"


class PixiInstallPasaTests(unittest.TestCase):
    def _write_fake_pasa_install(self, root):
        install_sh = root / "PASApipeline" / "scripts" / "install.sh"
        install_sh.parent.mkdir(parents=True)
        install_sh.write_text(
            textwrap.dedent(
                """\
                #!/usr/bin/env bash
                set -euo pipefail

                install_prefix=""
                while [ "$#" -gt 0 ]; do
                    case "$1" in
                        --install-prefix)
                            install_prefix="$2"
                            shift 2
                            ;;
                        *)
                            shift
                            ;;
                    esac
                done

                mkdir -p "${install_prefix}/bin" "${install_prefix}/src"
                for binary in pasa_rust slclust_rust cdbyank_rust faidx_rust; do
                    printf '#!/usr/bin/env bash\nexit 0\n' > "${install_prefix}/bin/${binary}"
                    chmod +x "${install_prefix}/bin/${binary}"
                done
                printf '#!/usr/bin/env perl\nprint qq(fake\\n);\n' > "${install_prefix}/src/Launch_PASA_pipeline.pl"
                chmod +x "${install_prefix}/src/Launch_PASA_pipeline.pl"
                ln -s ../bin "${install_prefix}/src/bin"
                """
            )
        )
        install_sh.chmod(install_sh.stat().st_mode | stat.S_IXUSR)
        return install_sh

    def _run_wrapper(self, cwd, conda_prefix):
        env = os.environ.copy()
        env["CONDA_PREFIX"] = str(conda_prefix)
        return subprocess.run(
            ["bash", str(PIX_INSTALL_PASA)],
            cwd=str(cwd),
            env=env,
            capture_output=True,
            text=True,
        )

    def test_local_install_uses_pasahome_layout(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            workdir = root / "work"
            conda_prefix = workdir / "env"
            workdir.mkdir()
            self._write_fake_pasa_install(root)

            result = self._run_wrapper(workdir, conda_prefix)

            self.assertEqual(
                result.returncode,
                0,
                msg=result.stdout + "\n" + result.stderr,
            )
            self.assertTrue(
                (conda_prefix / "opt" / "pasa" / "src" / "Launch_PASA_pipeline.pl").is_file()
            )
            self.assertTrue((conda_prefix / "opt" / "pasa" / "src" / "bin").exists())
            self.assertTrue((conda_prefix / "opt" / "pasa" / "bin" / "pasa_rust").is_file())

    def test_wrapper_is_idempotent_for_completed_install(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            workdir = root / "work"
            conda_prefix = workdir / "env"
            workdir.mkdir()
            pasa_repo = root / "PASApipeline"
            self._write_fake_pasa_install(root)

            first = self._run_wrapper(workdir, conda_prefix)
            self.assertEqual(first.returncode, 0, msg=first.stdout + "\n" + first.stderr)

            for child in list(pasa_repo.iterdir()):
                if child.is_dir():
                    shutil.rmtree(child)
                else:
                    child.unlink()

            second = self._run_wrapper(workdir, conda_prefix)
            self.assertEqual(second.returncode, 0, msg=second.stdout + "\n" + second.stderr)


if __name__ == "__main__":
    unittest.main()
