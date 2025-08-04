from distutils.core import setup

setup(
     name='leader_follower_KS',
     version='1.0',
     author='Juan Fernandez Gracia, Nils Kreuter',
     packages=['leader_follower_KS'],
     scripts=[],
     url='blank',
     license='LICENSE.txt',
     description='A package for detecting leader-follower dynamics in acoustic telemetry data.',
     install_requires=[
         "python >= 3.7.0",
         "datetime >= 5.1",
         "pandas >= 1.5.3",
         "matplotlib >= 3.7.1",
         "seaborn >= 0.12.2",
         "networkx >= 2.8.4",
         "scipy >= 1.11.4",
         "pytz >= 2022.7",
         "numpy >= 1.24.3"
     ],
     classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
