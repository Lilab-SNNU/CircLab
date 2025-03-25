import os
import sys
from os import path
script_path = path.abspath(__file__)
install_dir = path.dirname(os.path.dirname(script_path))
config_path = path.join(install_dir, 'config', 'config.yaml')