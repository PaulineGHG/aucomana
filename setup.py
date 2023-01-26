from setuptools import setup

setup(
    name='aucome_analysis',
    version='0.1.0',
    description='Python package furnishing library to analyse AuCoMe tool runs.',
    url='https://github.com/PaulineGHG/analysis_runs',
    author='Pauline HAMON-GIRAUD',
    license='GPL-3.0 license',
    packages=['aucome_analysis'],
    install_requires=['biopython>=1.79',
                      'matplotlib>=3.5.1',
                      'numpy>=1.22.3',
                      'pandas>=1.4.2',
                      'pyahocorasick>=1.4.4',
                      'rpy2>=3.5.1',
                      'venn>=0.1.3,'
                      ],

    classifiers=[
        'Development Status :: 1 - Beta',
        'Intended Audience :: Science/Research/Bio-informatics',
        'License :: OSI Approved',
        'Operating System :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        'Programming Language :: Python :: 3.7',
    ], 
  entry_points={
      'console_scripts': [
          'aucome_analysis = aucome_analysis.__main__:main']
)
