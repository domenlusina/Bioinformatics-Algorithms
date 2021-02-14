import gzip
import os
import zlib
from subprocess import Popen

import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform


def loadFasta(filename, verbose=0):
    """ Parses a classically formatted and possibly
        compressed FASTA file into a dictionary where the key
        for a sequence is the first part of its header without
        any white space; if verbose is nonzero then the identifiers
        together with lengths of the read sequences are printed"""
    if filename.endswith(".gz"):
        fp = gzip.open(filename, 'rt')
    else:
        fp = open(filename, 'r')
    # split at headers
    # data = fp.read().split('>')
    data = fp.read()
    data = data.split('>')
    fp.close()
    # ignore whatever appears before the 1st header
    data.pop(0)
    # prepare the dictionary
    D = {}
    for sequence in data:
        lines = sequence.split('\n')
        header = lines.pop(0).split()
        key = header[0]
        D[key] = ''.join(lines)
        if verbose:
            print("Sequence %s of length %d read" % (key, len(D[key])))
    return D


def draw_dendrogram(distance_matrix, txt_file_name, labels):
    """
    A function that draws a dendrogram for a given distance matrix and also saves distance matrix to a file.

    :param distance_matrix: distance matrix NxN of floats
    :param txt_file_name: file name that we wish to save distance table to
    :param labels: string array, name of a leaf in dendrogram
    :return:
    """
    distances = squareform(distance_matrix)
    linkage_matrix = linkage(distances, "single")

    np.savetxt(txt_file_name, distance_matrix, fmt='%1.3f')

    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('organisms')
    plt.ylabel('distance')
    dendrogram(
        linkage_matrix,
        labels=labels,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=8.,  # font size for the x axis labels
    )
    plt.show()


def generate_files(data, joined):
    """
    Generates files of DNA as a string and saves them. We also save concated strings.

    :param data: dictionary where items are represented as a string
    :param joined: dictionary of joined data items
    :return: 
    """
    print("Generating files of sequences.")
    for key in data:
        file = open("GenCompress\\" + key + ".dat", "w")
        file.write(data[key])
        file.close()

    for (key1, key2) in joined:
        file = open("GenCompress\\" + key1 + key2 + ".dat", "w")
        file.write(joined[(key1, key2)])
        file.close()


# reading data
data = loadFasta("Seq.fasta", 0)

# we sort organisms
keys = list(data.keys())
keys.sort()

# generating concated strings
joined = {}
n_org = len(keys)
for i in range(n_org):
    for j in range(i + 1, n_org):
        joined[(keys[i], keys[j])] = data[keys[i]] + data[keys[j]]

# change to True if we haven't generate .dat files yet
if False:
    # we save strings as .dat files
    generate_files(data, joined)
    # we run .bat file so we obtain output files
    p = Popen("GenCompress\\compress_all.bat")
    stdout, stderr = p.communicate()

GenCompress = True
MyCompress = True
if GenCompress:
    # obtaining length of compressed files from .out files
    compressed = {}
    for f in os.listdir('GenCompress\\out'):
        with open("GenCompress\\out\\" + f) as file:
            content = file.readlines()
            for line in content:
                if line.startswith(' The size of compressed file is'):
                    compressed[f[:-4]] = int(line.split()[-2])
                    break
            else:
                print(f)
        file.close()

    # calculating distance matrix
    distance_matrix = np.zeros((n_org, n_org))
    for i in range(n_org):
        for j in range(i + 1, n_org):
            d = 1 - (compressed[keys[i]] - compressed[keys[i] + "_" + keys[j]]) / compressed[keys[i] + keys[j]]
            distance_matrix[i, j] = d
            distance_matrix[j, i] = d

    # drawing a dendrogram
    draw_dendrogram(distance_matrix, "res_gen_compress.txt", keys)

if MyCompress:
    # compressing the data
    compressed = {}
    for key, item in data.items():
        compressed[key] = len(zlib.compress(bytes(item, "utf8")))

    for key, item in joined.items():
        compressed[key] = len(zlib.compress(bytes(item, "utf8")))

    # calculating the distance matrix
    distance_matrix = np.zeros((n_org, n_org))
    for i in range(n_org):
        for j in range(i + 1, n_org):
            d = 1 - (compressed[keys[i]] - (compressed[(keys[i], keys[j])] - compressed[keys[j]])) / compressed[
                (keys[i], keys[j])]
            distance_matrix[i, j] = d
            distance_matrix[j, i] = d
    # drawing a dendrogram
    draw_dendrogram(distance_matrix, "res_zlib_compress.txt", keys)
