from distutils.core import setup

setup(
    name=‘f3',
    version='0.1dev',
    packages=[‘photometry’,],
    license=‘MIT’,
    long_description=open('README.txt').read(),
    author=‘Ben Montet’,
    author_email=‘bmontet@uchicago.edu’,
    url=‘http://github.com/benmontet/f3’,
    description=‘Photometry for Kepler Full Frame Images’,
    install_requires=[
          ‘mahotas’,
      ],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python 2.7”,
    ],
)
