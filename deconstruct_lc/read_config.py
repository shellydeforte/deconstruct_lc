import configparser
import os


def read_config():
    config = configparser.ConfigParser()
    cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                       'config.cfg'))
    config.read_file(open(cfg_fp, 'r'))
    return config


def read_test_config():
    config = configparser.ConfigParser()
    cfg_fp = os.path.join(os.path.join(os.path.dirname(__file__), '..',
                                       'config_test.cfg'))
    config.read_file(open(cfg_fp, 'r'))
    return config


def main():
    pass


if __name__ == '__main__':
    main()