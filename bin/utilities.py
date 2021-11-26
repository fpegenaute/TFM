import os


def get_filename_ext(filepath):
    """
    Given a file path, split it by dots and get the extension as the last element
    Also, take the basename in the path and split it by the point, get the first
    element as the name of the file
    """

    extension = filepath.split(".")[-1]
    filename = os.path.basename(filepath).split('.')[0]
    
    return filename, extension


