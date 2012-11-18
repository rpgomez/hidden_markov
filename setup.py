from numpy.distutils.core import setup

from numpy.distutils.misc_util import Configuration

def configuration(parent_package='',top_path=None):
    config = Configuration('hiddenmarkov',parent_package,top_path)
    config.add_extension('hiddenmarkov', ['hidden_markov.f90'])
    return config

setup(configuration=configuration)
