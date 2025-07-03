from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="sprint-satellite",
    version="1.0.0",
    author="Sprint Satellite Project",
    author_email="your-email@domain.com",
    description="Simulation of high-power, short-duration CubeSat missions using electrodynamic tethers",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/sprint-satellite",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    python_requires=">=3.7",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "sprint-sat-sim=sprint_sat_sim:main",
        ],
    },
    keywords="satellite, cubesat, mhd, electrodynamic, tether, orbital, simulation",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/sprint-satellite/issues",
        "Source": "https://github.com/yourusername/sprint-satellite",
        "Documentation": "https://github.com/yourusername/sprint-satellite#readme",
    },
) 