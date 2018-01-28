
import os
import numpy as np

def walk_up_folder(path, depth=1):
    """ Walk up a specific number of sub-directories """
    _cur_depth = 0
    while _cur_depth < depth:
        path = os.path.dirname(path)
        _cur_depth += 1
    return path

def remove_fields(nparray, names):
    """ Remove fields field from a structured numpy array """

    fields = list(nparray.dtype.names)
    for name in names:
        if name in fields:
            fields.remove(name)
    return nparray[fields]


def find_files(path,searchstr):
    """ Recursive file search with a given search string """

    res = []
    for root, dirs, files in os.walk(path):
        for f in files:
            if f.find(searchstr) != -1:
                res.append(os.path.join(root, f))

    if len(res) == 0:
        print 'No files found which contain: "' + param + '".'
        return None
    elif len(res) == 1:
        return res[0]
    else:
        return np.array(res)

