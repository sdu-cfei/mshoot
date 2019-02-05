from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()  # TODO: Images might not work on PyPi

setup(name='mshoot',
      version='0.0.2dev',
      description='Multiple shooting dynamic optimization',
      long_description=long_description,
      url='https://github.com/sdu-cfei/mshoot',
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
        'pandas',
        'scipy',
        'matplotlib',
        'scikit-learn',
        'pyfmi'
      ],
      classifiers = [
          'Programming Language :: Python :: 3',
          'Topic :: Scientific/Engineering',
          'License :: OSI Approved :: BSD License'
      ]
)