import sys


def test_unicode(f):
    """
    attempt to open input file as unicode
    if file fails to be read using readline method, assume .gz
    compression
    """
    try:
        with open(f) as f:
            f.readline()
        compressed = False
    except UnicodeDecodeError:
        compressed = True
    except IsADirectoryError:
        return None

    return compressed


if __name__ == '__main__':
    f = sys.argv[1]
    result = test_unicode(f)
    if result:
        print(f'{f} is gzipped')
    else:
        print(f'{f} is not gzipped')

