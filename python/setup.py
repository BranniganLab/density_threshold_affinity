from setuptools import setup

setup(
    name='DTA',
    version='1.0',
    description='Tools for Doing Density-Threshold Affinity analysis',
    url='https://github.com/BranniganLab/density_threshold_affinity',
    author='Brannigan Lab',
    author_email='grace.brannigan@rutgers.edu',
    packages=['DTA'],
    install_requires=['numpy>=2.0.0','pandas>=2.2.2', 'matplotlib>=3.9.1', 'scipy>=1.14.0'],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5'
    ],
)
