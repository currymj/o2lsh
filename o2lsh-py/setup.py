import sys

from setuptools import setup
from setuptools.command.test import test as TestCommand

try:
    from setuptools_rust import RustExtension
except ImportError:
    import subprocess
    errno = subprocess.call([sys.executable, '-m', 'pip', 'install', 'setuptools-rust'])
    if errno:
        print("Please install setuptools-rust package")
        raise SystemExit(errno)
    else:
        from setuptools_rust import RustExtension


class PyTest(TestCommand):
    user_options = []

    def run(self):
        self.run_command("test_rust")

        import subprocess
        import sys
        errno = subprocess.call([sys.executable, '-m', 'pytest', 'tests'])
        raise SystemExit(errno)


setup_requires = ['setuptools-rust>=0.6.1']
install_requires = []
tests_require = install_requires + ['pytest', 'pytest-benchmark']

setup(
    name='o2lshpy',
    version='0.1.0',
    packages=['o2lsh'],
    rust_extensions=[RustExtension('o2lsh._o2lsh', 'Cargo.toml')],
    install_requires=install_requires,
    tests_require=tests_require,
    setup_requires=setup_requires,
    include_package_data=True,
    zip_safe=False
)
