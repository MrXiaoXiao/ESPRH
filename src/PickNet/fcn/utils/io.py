import os
import yaml
from time import strftime, localtime
from termcolor import colored


class IO():

    def __init__(self, log_dir=None):

        self.log_dir = log_dir

    def read_yaml_file(self, config_file):

        pfile = open(config_file)
        d = yaml.load(pfile)
        pfile.close()

        return d

    def print_info(self, info_string, quite=False):

        info = '[{0}][INFO] {1}'.format(self.get_local_time(), info_string)
        print(colored(info, 'green'))

    def print_warning(self, warning_string):

        warning = '[{0}][WARNING] {1}'.format(self.get_local_time(), warning_string)

        print(colored(warning, 'blue'))

    def print_error(self, error_string):

        error = '[{0}][ERROR] {1}'.format(self.get_local_time(), error_string)

        print(colored(error, 'red'))

    def get_local_time(self):

        return strftime("%d %b %Y %Hh%Mm%Ss", localtime())
