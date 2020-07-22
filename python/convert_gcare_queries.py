import sys
from os import listdir
from os.path import isfile, join

query_dir = sys.argv[1]
true_dir = sys.argv[2]
query_files = [f for f in listdir(query_dir) if isfile(join(query_dir, f))]
true_files = [f for f in listdir(true_dir) if isfile(join(true_dir, f))]

true_cards = {}
queries = []

def main():
    for true_file in true_files:
        with open(true_dir + '/' + true_file, 'r') as f:
            for line in f:
                card = line.strip()
                if true_file not in true_cards:
                    true_cards[true_file] = card
    

    for query_file in query_files:
        with open(query_dir + '/' + query_file, 'r') as f:
            query = []
            for line in f:
                if line.startswith('e'):
                    edge = line.strip().split()
                    query.append((edge[1], edge[3], edge[2]))
            queries.append((query, true_cards[query_file]))
    
    
    with open(sys.argv[3], 'w') as f:
        for query, card in queries:
            v_list, label_seq = format_query(query)
            f.write(v_list + ',' + label_seq + ',' + card + '\n')


def format_query(query):
    srcdest_pair = []
    label = {}
    for edge in query:
        srcdest_pair.append((edge[0], edge[2]))
        label[(edge[0], edge[2])] = edge[1]

    srcdest_pair.sort()
    v_list = []
    label_seq = []
    for srcdest in srcdest_pair:
        v_list.append(srcdest[0] + '-' + srcdest[1])
        label_seq.append(label[srcdest])

    return ';'.join(v_list), '->'.join(label_seq)



if __name__ == '__main__':
    main()

