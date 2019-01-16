from setuptools import setup, find_packages

setup(
    name="crosssection",
    version="0.1",
    description="A totally irrelevant tool with no practical use case.",
    author="Kyle Poland",
    packages=find_packages(),
    install_requires=["mpmath", "numpy", "scipy", "matplotlib"],
    zip_safe=False,
    entry_points={"console_scripts": ["crosssection=crosssection.main:main"]},
)