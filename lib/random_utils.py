from random import sample

def partition_list(list, k):
    """
    Partitions the list into k partitions returning a list of tuples with the start and end indices of each partition
    """
    if k > len(list):
        raise Exception("K must be smaller than the length of the list")
    base_points_per_partition = len(list) // k
    points_for_last_partition = base_points_per_partition + (len(list) - (base_points_per_partition * k))

    def points_per_partition(i):
        if i == k - 1:
            return points_for_last_partition
        return base_points_per_partition

    partitions = [(i * base_points_per_partition, (i * base_points_per_partition) + points_per_partition(i)) for i in range(k)]
    return partitions

def random_n_elements_across_k_partitions(list, n, k):
    """
    Returns n random elements from the list, with the constraint that
    the elements are well spread out across k partititions of the list

    Note: You may not get N points if N is not divisible by k
    """
    if k > n:
        raise Exception("N must be larger than k")
    if n > len(list):
        raise Exception("N must be smaller than the length of the list")

    base_points_per_partition = n // k
    points_for_last_partition = base_points_per_partition + (n - (base_points_per_partition * k))

    def points_per_partition(i):
        if i == k - 1:
            return points_for_last_partition
        return base_points_per_partition

    partitions = partition_list(list, k) 
    
    random_indices = [sample(range(partitions[i][0], partitions[i][1]), points_per_partition(i)) for i in range(k)]
    flattened = [item for sublist in random_indices for item in sublist]
    return flattened 


def test_partition_list():
    list = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    k = 4
    partitions = partition_list(list, k)
    assert(len(partitions) == k)
    assert(partitions[0] == (0, 2))
    assert(partitions[1] == (2, 4))
    assert(partitions[2] == (4, 6))
    assert(partitions[3] == (6, 9))

    list = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    k = 2
    partitions = partition_list(list, k)
    assert(len(partitions) == k)
    assert(partitions[0] == (0, 4))
    assert(partitions[1] == (4, 9))

    list = [1, 2, 3]
    k = 1
    partitions = partition_list(list, k)
    assert(len(partitions) == k)
    assert(partitions[0] == (0, 3))


def test_random_n_elements_from_list():
    # You get exactly N points 
    list = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    n = 6
    k = 4
    indices = random_n_elements_across_k_partitions(list, n, k)
    assert(len(indices) == n)


    # Works for cases when number of points (n) is not divisible by number of partitions (k)
    list = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    n = 7
    k = 2
    
    indices = random_n_elements_across_k_partitions(list, n, k)
    assert(len(indices) == n)
    # No duplicate points
    
