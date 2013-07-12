from distutils.core import setup

setup(name='trep_collisions',
      version='0.1',
      description='solves collisions in trep',
      author='Elliot Hevel',
      author_email='elliothevel2013@u.northwestern.edu',
      package_dir={'trep_collisions': 'src'},
      packages=['trep_collisions', 'trep_collisions.surfaces'])
