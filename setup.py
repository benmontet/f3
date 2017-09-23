from setuptools import setup
import f3

setup(
    name='f3',
    version='1.0.4',
    license='MIT',
    long_description=open('README.rst').read(),
    author='Ben Montet',
    author_email='bmontet@uchicago.edu',
    packages=[
        'f3',
    ],
    include_package_data=True,
    url='http://github.com/benmontet/f3',
    description='Photometry for Kepler Full Frame Images',
    package_data={'': ['README.rst', 'LICENSE']},
    install_requires=[
          'mahotas',
      ],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
    ],
)
