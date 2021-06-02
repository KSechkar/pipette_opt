import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ppopt",
    version="0.1.0",
    author="Kirill Sechkar",
    author_email="kirill.sechkar18@imperial.ac.uk",
    description="Saving pipette tips when distributing DNA parts in automated DNA assembly",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/KSechkar/pipette_opt",
    project_urls={
        "Bug Tracker": "https://github.com/KSechkar/pipette_opt/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.8",
    install_requires=['numpy>=1.17.2','gurobipy>=9.1.1', 'opentrons==3.21.1']
)