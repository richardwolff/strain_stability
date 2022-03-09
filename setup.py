from setuptools import setup

setup(
   name='strainstability',
   version='0.1.0',
   author='Richard Wolff',
   author_email='rwolff@g.ucla.edu',
   packages=['package_name', 'package_name.test'],
   description='Scripts for analyzing the stability of strains',
   long_description=open('README.txt').read(),
   install_requires=[
       "Django >= 1.1.1",
       "pytest",
   ],
)