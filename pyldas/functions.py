
import os

def remove_fields(nparray, names):
    fields = list(nparray.dtype.names)
    for name in names:
        if name in fields:
            fields.remove(name)
    return nparray[fields]


def find_files(path,name):
    res = []
    for root, dirs, files in os.walk(path):
        for f in files:
            if f.find(name) != -1:
                res.append(os.path.join(root, f))

    if len(res) == 0:
        return None
    elif len(res) == 1:
        return res[0]
    else:
        return res
