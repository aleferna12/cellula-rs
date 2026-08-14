import numpy as np


def bitarray_to_dec(a):
    return int("".join(a.astype(int).astype(str)), 2)


def dec_to_bitarray(num, width=None):
    return np.array(list(np.binary_repr(num, width)), dtype=int)


def calculate_complimentarity(key, lock):
    return np.bitwise_count(key ^ lock)


def cell_contact_energy(k1: int, l1: int, k2: int, l2: int, cell_energy, gene_energy, genome_length):
    comp = calculate_complimentarity(k1, l2) + calculate_complimentarity(k2, l1)
    return 2 * (cell_energy + gene_energy * (0.5 - comp / (2 * genome_length)))


def calculate_gamma(med_energy, contact_energy):
    return med_energy - contact_energy / 2


def min_width(bitdec):
    bitdec = np.where(bitdec == 0, 1, bitdec)
    return np.ceil(np.log2(bitdec + 1)).astype(int)


def bitstr_to_bitdec(keydec, lockdec, keylock_width=None):
    if keylock_width is None:
        keylock_width = min_width(max(keydec, lockdec))
    binary_a = np.binary_repr(keydec, width=keylock_width)  # Convert decimal to binary string
    binary_b = np.binary_repr(lockdec, width=keylock_width)
    concatenated_binary = binary_a + binary_b
    return int(concatenated_binary, 2)


def hamming_distance(bitstr1, bitstr2, sep="-"):
    k1, l1 = np.fromstring(bitstr1, sep=sep, dtype=int)
    k2, l2 = np.fromstring(bitstr2, sep=sep, dtype=int)
    width = np.max(min_width([k1, l1, k2, l2]))
    bitdec1 = bitstr_to_bitdec(k1, l1, width)
    bitdec2 = bitstr_to_bitdec(k2, l2, width)
    # Perform bitwise XOR operation
    xor_result = bitdec1 ^ bitdec2

    # Count the number of set bits (1s)
    hamming_dist = 0
    while xor_result:
        hamming_dist += xor_result & 1
        xor_result >>= 1

    return hamming_dist
