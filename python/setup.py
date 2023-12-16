from setuptools import setup

setup(
    name='DTA',
    version='0.2',
    description='Tools for Doing Density Threshold Affinity Analysis',
    url='https://github.com/BranniganLab/Polar_Binning_DeltaG',
    author='Brannigan Lab',
    author_email='grace.brannigan@rutgers.edu',
    packages=['DTA'],
    install_requires=['numpy>=1.22.0','pandas>=1.4.0', 'matplotlib>=3.5.0'],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5'
    ],
)
