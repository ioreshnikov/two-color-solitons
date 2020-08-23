from setuptools import setup


setup(
    name="common",
    version="0.1",
    description=(
        "The simulation and the plotting codes for the paper "
        "on two-color solitons."),
    author="Ivan Oreshnikov",
    author_email="oreshnikov.ivan@gmail.com",
    python_requires=">=3.5",
    install_requires=[
        "matplotlib==3.3.0",
        "numpy==1.19.1",
        "scipy==1.5.2",
        "sympy==1.6.1"
    ],
    scripts=[]
)
