'''
Find proto-restriction sites in an arbitrary DNA sequence, given a list of sites

arguments: Master RE list, RE search list, DNA file

written for Python 3.5.2
'''
import os, sys


# creats a dictionary of RE sites from the file passed
def load_sites(file, sites, cut_site=0):
    site_dict = {}
    with open(file, 'r', encoding='utf-8') as a:
        temp = a.read()
    
    temp = temp[1:] #remove leading utf char
    temp = temp.split('\n')
    
    # The next three statements remove duplicates (list->set), and
    # removes empty elements (from blank lines at the EOF.
    temp = set(temp)
    temp = list(temp)
    try:
        temp.remove('')
    except:
        pass
    
    for i in range(0, len(temp)):
        temp[i] = temp[i].split(',')
        if cut_site==0:
            temp[i][1] = temp[i][1].replace('^','')
    
    temp = dict(temp)
    
    with open(sites, 'r', encoding='utf-8') as b:
        site_list = b.read()
    
    site_list = site_list[1:]
    site_list = site_list.split('\n')
    site_list = list(set(site_list))
    try:
        site_list.remove('')
    except:
        pass
        
    for i in range(0, len(site_list)):
        site_dict[site_list[i]] = temp[site_list[i]]
    
    return site_dict
    

# Reads the plasmid sequence 
def load_DNA(file):
    with open(file, 'r', encoding='utf-8') as a:
        seq = a.read()
    
    #seq = seq[1:] # remove leading utf char
    return seq


def create_search_matrix(site):
    matrix = blank_matrix()
    
    for i in range(0, len(site)):
        matrix[get_index(site[i])][i] = 1
    
    return matrix


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


def blank_matrix():
    matrix = [[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]]
    return matrix


def score_window(window, matrix):
    score = 0
    for i in range(0, len(window)):
        a = window[i]
        b = get_index(a)
        
        score += matrix[b][i]
    
    return score


if __name__ == '__main__':
    master, search, dna_file = sys.argv[1], sys.argv[2], sys.argv[3]

    site_dict = load_sites(master, search)
    site_keys = [s for s in site_dict.keys()]
    dna = load_DNA(dna_file)
    
    l = len(dna)
    dna += dna[0:6]
    threshold = 5
    stuff = 'Enzyme,Pos,Score,Seq\n'

    for j in range(0, len(site_dict)):
        matrix = create_search_matrix(site_dict[site_keys[j]])
        for i in range(0,l):
            window = dna[i:i+6]
            a = score_window(window, matrix)
            if a >= threshold:
                #print('{3} {0}({1}): {2}'.format(i+1, a, window, site_keys[j]))
                stuff += ','.join([site_keys[j], str(i+1), str(a), window])
                stuff += '\n'
                
    filename = os.path.splitext(dna_file)[0] + " site list.csv"
    with open(filename, 'w', encoding='utf-8') as csv:
        csv.write(stuff)
        










