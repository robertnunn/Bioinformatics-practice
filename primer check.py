'''
Given a list of sites and a sequence file to find them in, find all sites that have above a
user-defined threshold percentage identity

arguments: DNA file, search list, threshold %

written for Python 3
'''
import os, sys


# creates a list of sites from the file passed
def load_sites(sites):
    with open(sites, 'r', encoding='utf-8') as a:
        temp = a.read()

    temp = temp[1:] # remove leading utf char
    temp = temp.split('\n')

    # The next three statements remove duplicates (list->set), and
    # removes empty elements (from blank lines at the EOF.
    temp = set(temp)
    temp = list(temp)
    try:
        temp.remove('')
    except:
        pass

    return temp


# Reads the plasmid sequence
def load_DNA(file):
    with open(file, 'r', encoding='utf-8') as a:
        seq = a.read()

    # seq = seq[1:] # remove leading utf char
    seq = seq.replace('\n', '')
    return seq.upper()


def create_search_matrix(site):
    score = dict()
    # single letter nucleotide codes
    # [A,C,G,T]
    score["A"] = [1,0,0,0]
    score["B"] = [0,1,1,1]
    score["C"] = [0,1,0,0]
    score["D"] = [1,0,1,1]
    score["G"] = [0,0,1,0]
    score["H"] = [1,1,0,1]
    score["K"] = [0,0,1,1]
    score["M"] = [1,1,0,0]
    score["N"] = [1,1,1,1]
    score["R"] = [1,0,1,0]
    score["S"] = [0,1,1,0]
    score["T"] = [0,0,0,1]
    score["U"] = [0,0,0,1]
    score["V"] = [1,1,1,0]
    score["W"] = [1,0,0,1]
    score["X"] = [1,1,1,1]
    score["Y"] = [0,1,0,1]

    search_matrix = []

    for i in site:
        search_matrix.append(score[i])

    return search_matrix


def get_index(char):
    if char == 'A':
        return 0
    elif char == 'C':
        return 1
    elif char == 'G':
        return 2
    elif char == 'T':
        return 3
    else:
        return -1


def score_window(window, matrix):
    score = 0
    for i in range(len(window)):
        a = window[i]
        b = get_index(a)

        score += matrix[i][b]

    return score


def make_alignment(window, site):
    align = ''
    for i in range(len(window)):
        if window[i] != site[i]:
            align += ' '
        else:
            align += "|"

    return align


if __name__ == '__main__':
    dna_file, search, threshold = sys.argv[1], sys.argv[2], sys.argv[3]

    # dna_file = "E coli rrsA operon.txt"
    # search = "probes.txt"
    # threshold = 55
    threshold = float(threshold)
    threshold /= 100.0

    site_list = load_sites(search)
    dna = load_DNA(dna_file)

    max_len = 0
    for i in site_list:
        if len(i) > max_len:
            max_len = len(i)

    l = len(dna)
    dna += dna[0:max_len]
    stuff = 'Query,Pos,Score,% Identity,Subject,Alignment,Threshold=' + str(threshold) + '\n'

    for j in site_list:
        matrix = create_search_matrix(j)
        for i in range(l):
            window = dna[i:i+len(j)]
            a = score_window(window, matrix)
            b = a / len(j)
            if b >= threshold:
                # print('{3} {0}({1}): {2}'.format(i+1, a, window, j))
                c = j + "\n,,,,," + make_alignment(window, j) + "\n,,,,," + window
                stuff += ','.join([j, str(i+1), str(a), str(round(b, 3)), window, c])
                stuff += '\n'

    filename = os.path.splitext(dna_file)[0] + " site list.csv"
    with open(filename, 'w', encoding='utf-8') as csv:
        csv.write(stuff)
