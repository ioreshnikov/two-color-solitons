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
        "humanize",
        "matplotlib==3.3",
        "numpy",
        "scipy",
        "sympy"
    ],
    scripts=[]
)
