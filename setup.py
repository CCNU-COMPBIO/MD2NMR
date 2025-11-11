from setuptools import setup, find_packages

setup(
    name="MD2NMR",
    version="0.0.6",
    author="Tiejun Wei and Houfang Zhang",
    description="Tools for calculating NMR relaxation observables (R1/R2/NOE/T1/T2/Tau_c) directly from MD trajectories. Initially written for calculations regarding nucleosome simulations but can be extended for other proteins/complexes.",
    url="https://github.com/HoufangZhang/MD2NMR",
    packages=find_packages(),
    install_requires=["numpy == 2.2.6",
                      "pandas == 2.3.3",
                      "scikit-learn == 1.7.2",
                      "scipy == 1.15.3",
                      "MDAnalysis == 2.9.0",
                      "matplotlib == 3.4.3",
                      'setuptools == 65.5.0'],
)
