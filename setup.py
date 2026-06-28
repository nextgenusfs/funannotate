#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Note: To use the 'upload' functionality of this file, you must:
#   $ pip install twine

import io
import os
import sys
from shutil import rmtree

from setuptools import find_packages, setup, Command
from setuptools.command.build_py import build_py as _build_py

# Package meta-data.
NAME = "funannotate"
DESCRIPTION = "funannotate: eukaryotic genome annotation pipeline"
URL = "https://github.com/nextgenusfs/funannotate"
EMAIL = "nextgenusfs@gmail.com"
AUTHOR = "Jon Palmer"
REQUIRES_PYTHON = ">=3.6.0, <3.12"
VERSION = None

# What packages are required for this module to be executed?
REQUIRED = [
    "biopython<1.80",
    "goatools",
    "seaborn",
    "psutil",
    "pandas",
    "matplotlib",
    "natsort",
    "numpy",
    "requests",
    "scikit-learn",
    "scipy",
    "distro",
    "packaging",
]

# What packages are optional?
EXTRAS = {
    # 'fancy feature': ['django'],
}

# The rest you shouldn't have to touch too much :)
# ------------------------------------------------
# Except, perhaps the License and Trove Classifiers!
# If you do change the License, remember to change the Trove Classifier for that!

here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
try:
    with io.open(os.path.join(here, "README.md"), encoding="utf-8") as f:
        long_description = "\n" + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

# Load the package's __version__.py module as a dictionary.
# __file__ must be injected so _git_version() can locate the .git directory.
about = {"__file__": os.path.join(here, NAME, "__version__.py")}
if not VERSION:
    with open(os.path.join(here, NAME, "__version__.py")) as f:
        exec(f.read(), about)
else:
    about["__version__"] = VERSION

# Bake the resolved version into _version.txt so eggs/wheels without .git
# can still report the full version string at runtime.
with open(os.path.join(here, NAME, "_version.txt"), "w") as _vf:
    _vf.write(about["__version__"] + "\n")


class UploadCommand(Command):
    """Support setup.py upload."""

    description = "Build and publish the package."
    user_options = []

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print(("\033[1m{0}\033[0m".format(s)))

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.status("Removing previous builds…")
            rmtree(os.path.join(here, "dist"))
        except OSError:
            pass

        self.status("Building Source and Wheel (universal) distribution…")
        os.system("{0} setup.py sdist bdist_wheel --universal".format(sys.executable))

        self.status("Uploading the package to PyPI via Twine…")
        os.system("twine upload dist/*")

        self.status("Pushing git tags…")
        os.system("git tag v{0}".format(about["__version__"]))
        os.system("git push --tags")

        sys.exit()


class BuildPyCommand(_build_py):
    """build_py that preserves the execute bit on bundled helper scripts.

    The scripts in ``funannotate/aux_scripts`` are shipped as package data and
    several of them are launched directly via their shebang at runtime (the
    Perl/shell helpers, phobius-multiproc.py, etc.). When funannotate is
    installed from a wheel/sdist the copied data files do not reliably keep
    their execute bit, leaving those scripts non-executable after install and
    breaking the pipeline (issue #1129). Re-assert the execute bits on the
    build output so the wheel -- and therefore the installation -- keep them.
    """

    _EXECUTABLE_DIRS = (os.path.join("funannotate", "aux_scripts"),)

    def run(self):
        super().run()
        for reldir in self._EXECUTABLE_DIRS:
            target = os.path.join(self.build_lib, reldir)
            if not os.path.isdir(target):
                continue
            for name in os.listdir(target):
                path = os.path.join(target, name)
                if os.path.isfile(path):
                    mode = os.stat(path).st_mode
                    # mirror the read bits into execute bits (u/g/o), e.g. 0644 -> 0755
                    os.chmod(path, mode | ((mode & 0o444) >> 2))


# Where the magic happens:
setup(
    name=NAME,
    version=about["__version__"],
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(exclude=("tests",)),
    package_data={"funannotate": ["_version.txt"]},
    entry_points={
        "console_scripts": ["funannotate=funannotate.funannotate:main"],
    },
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    include_package_data=True,
    license="BSD-2",
    # scripts=['scripts/funannotate'],
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python",
        "Operating System :: Unix",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    cmdclass={
        "upload": UploadCommand,
        "build_py": BuildPyCommand,
    },
)
