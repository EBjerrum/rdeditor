from setuptools import setup

setup(name='rdeditor',
      version='0.1',
      description='An RDKit based molecule editor using PySide',
      long_description="""The editor is not meant as a drawing program, rather a simple and easy to expand editor""",
      keywords="RDKit molecule editor pyside",
      url='http://github.com/ebjerrum/rdeditor',
      author='Esben Jannik Bjerrum',
      author_email='esbenjannik@rocketmail.com',
      license='LGPL',
      packages=['rdeditor'],
      package_data={'rdeditor': ['pixmaps/*']},
      install_requires=[
	'PySide2',
	'numpy'
      ],
      entry_points = {
              'console_scripts': [
                  'rdEditor = rdeditor.rdEditor:launch',
              ], },
      zip_safe=False)

