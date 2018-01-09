import itertools


def ints_to_ranges(int_list):
    """
    Given a list of integers (usually corresponding to index locations in a
    list), return a generator that provides a list of ranges for the indexes.
    example: [1,2,3,4,6,7,8] yeilds (1, 4), (6, 8)
    """
    for a, b in itertools.groupby(enumerate(int_list), lambda x: x[0] - x[1]):
        b = list(b)
        yield b[0][1], b[-1][1]


def sort_list_by_second_list(sort_by, dep_sort, flag=True):
    """
    This expression is broken down for learning.
    lambda defines a function. Given the variable pair, return the 0 indexed value
    zip combines two lists into a list of tuples
    The sorted method will accept a list of tuples, and sort by a specific
    value in each tuple. This is specified by the key.
    """
    key_fun = lambda pair: pair[0]
    list_tups = zip(sort_by, dep_sort)
    sorted_tups = sorted(list_tups, reverse=flag, key=key_fun)
    sorted_by, dep_sorted = zip(*sorted_tups)
    return list(sorted_by), list(dep_sorted)


def main():
    pass


if __name__ == '__main__':
    main()