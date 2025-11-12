import sys

def reducer():
    _map = dict()                                   # it's a map structure, which is efficient to do the searching and insertion.
    for line in sys.stdin:
        word, count = line.split("\t")              # since the {key, value} pair is separated by a '\t'
                                                    # Here is to use '\t' as the position to determine where to cut

        count = int(count)                          # do the type transformation to integer
        word = word.lower()                         # make sure to parse lower case

        if word in _map:
            _map[word] += 1                         # if the {key, value} pair exists, increase the value by 1
        else:                                       # otherwise, insert the {key, value} pair
            _map[word] = 1

    for key, val in _map.items():                   # trace all elements in the map
        print(key, val)


if __name__ == "__main__":
    reducer()
