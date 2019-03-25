
import os
import logging

import numpy as np

logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)

def walk_up_folder(path, depth=1):
    """ Walk up a specific number of sub-directories """
    _cur_depth = 0
    while _cur_depth < depth:
        path = os.path.dirname(path)
        _cur_depth += 1
    return path

def find_files(path,searchstr):
    """ Recursive file search with a given search string """

    res = []
    for root, dirs, files in os.walk(path):
        for f in files:
            if f.find(searchstr) != -1:
                res.append(os.path.join(root, f))

    if len(res) == 0:
        logging.warning('No files found which contain: "' + searchstr + '".')
        return None
    elif len(res) == 1:
        return res[0]
    else:
        return np.array(res)
