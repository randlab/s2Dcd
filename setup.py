import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="s2Dcd", # Replace with your own username
    version="0.1.0",
    author="Alessandro Comunian",
    author_email="alessandro DOT comunian AT gmail DOT com",
    description="Sequential 2D MPS simulations with conditioning data",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/randlab/s2Dcd",
    packages=setuptools.find_packages(where="s2Dcd"),
#    package_dir={"", "s2Dcd"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU GPL License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    include_package_data = True,
    package_data = {
        'in_text' : ["deesse_file_in_text.json"]
    },
    install_requires=[
        "numpy>=1.18.1",
        "pandas>=0.25.3"
        ],
    py_modules = ["deesse", "utili", "grid", "ext", "gslibnumpy"],
)
