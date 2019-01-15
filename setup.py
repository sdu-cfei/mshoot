from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='mshoot',
      version='0.0.2dev',
      description='Multiple shooting dynamic optimization',
      long_description=long_description,
      url='https://github.com/sdu-cfei/mshoot',  # TODO: Not there yet...
      keywords='multiple shooting dynamic optimization',
      author='Krzysztof Arendt',
      author_email='krzysztof.arendt@gmail.com',
      license='BSD',
      platforms=['Windows', 'Linux'],
      packages=[
          'mshoot',
          'test'
      ],
      include_package_data=True,
      install_requires=[
          'numpy',
          'matplotlib',
          'scipy',
          'pandas',
          'dask',
      ],
      classifiers = [
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering',
          'License :: OSI Approved :: BSD License'
      ]
)