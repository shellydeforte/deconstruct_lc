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


def sort_lists(zip_list, flag=True):
    """
    Accepts a list of tuples created by zip(list1, list2, list3...)
    Will sort by the first item in the tuple
    lambda defines a function.
    Given the variable pair, return the 0 indexed value
    unpack either by for var1, var2, var3 in sorted_tups:
    or list1, list2, list3 = zip(*sorted_tups)
    """
    key_fun = lambda pair: pair[0]
    sorted_tups = sorted(zip_list, reverse=flag, key=key_fun)
    return sorted_tups


def demonstrate_sort():
    a = [2, 1, 3]
    b = ['cat', 'the', 'meowed']
    c = [6, 5, 7]
    list_tups = zip(a, b, c)
    key_fun = lambda pair: pair[0]
    sorted_tups = sorted(list_tups, reverse=True, key=key_fun)
    return sorted_tups


def main():
    demonstrate_sort()


if __name__ == '__main__':
    main()