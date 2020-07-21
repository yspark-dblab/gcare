import sys
import re
# write a function to decompose all the 2-paths from a given query

# for every pair, depending on dir combination, call corresponding function
 # ->->
 # <-->
 # -><-
# opt: if dir comb and labels are the same, see if we have computed it, then just re-use the cardinality

# have a graph structure of (label -> (src, dst)) so that we can iterate thru the edges, and compute degrees
 # also have graph structures using src/dst as key (to extend)

template_id = sys.argv[3]
query_len = 0
label2srcdest = {}
src2label2dst = {}
dst2label2src = {}
queries = []
cat2 = {}
cat1 = {}
pcat_max_deg = {}

def main():
    global query_len
    read_graph(sys.argv[1])
    read_queries(sys.argv[2])
    for query in queries:
        query_len = len(query)
        decoms = decompose(query)
        compute_cat_and_max_deg(decoms)
        compute_cat1(query)
    persist(sys.argv[4], sys.argv[5], sys.argv[6])

def compute_cat_and_max_deg(decoms):
    for decom in decoms:
        first = decom[0]
        sec = decom[1]

        v_list = first[0] + '-' + first[2] + ';' + sec[0] + '-' + sec[2]
        label_seq = first[1] + '->' + sec[1]
        if (v_list + ',' + label_seq) in cat2:
            continue

        card = 0
        max_deg = -1
        for f_src, f_dst in label2srcdest[first[1]]:
            if first[0] == sec[0]:
                if sec[1] in src2label2dst[f_src]:
                    deg = len(src2label2dst[f_src][sec[1]])
                    card += deg
                    max_deg = max(max_deg, deg)
            elif first[0] == sec[2]:
                if f_src in dst2label2src and sec[1] in dst2label2src[f_src]:
                    deg = len(dst2label2src[f_src][sec[1]])
                    card += deg
                    max_deg = max(max_deg, deg)
            elif first[2] == sec[0]:
                if f_dst in src2label2dst and sec[1] in src2label2dst[f_dst]:
                    deg = len(src2label2dst[f_dst][sec[1]])
                    card += deg
                    max_deg = max(max_deg, deg)
            elif first[2] == sec[2]:
                if sec[1] in dst2label2src[f_dst]:
                    deg = len(dst2label2src[f_dst][sec[1]])
                    card += deg
                    max_deg = max(max_deg, deg)

        cat2[v_list + ',' + label_seq] = card
        pcat_max_deg[first[0] + '-' + first[2] + ',' + first[1] + ',' + sec[0] + '-' + sec[2] + ',' + sec[1]] = max_deg

        card = 0
        max_deg = -1
        for s_src, s_dst in label2srcdest[sec[1]]:
            if sec[0] == first[0]:
                if first[1] in src2label2dst[s_src]:
                    deg = len(src2label2dst[s_src][first[1]])
                    card += deg
                    max_deg = max(max_deg, deg)
            elif sec[0] == first[2]:
                if s_src in dst2label2src and first[1] in dst2label2src[s_src]:
                    deg = len(dst2label2src[s_src][first[1]])
                    card += deg
                    max_deg = max(max_deg, deg)
            elif sec[2] == first[0]:
                if s_dst in src2label2dst and first[1] in src2label2dst[s_dst]:
                    deg = len(src2label2dst[s_dst][first[1]])
                    card += deg
                    max_deg = max(max_deg, deg)
            elif sec[2] == first[2]:
                if first[1] in dst2label2src[s_dst]:
                    deg = len(dst2label2src[s_dst][first[1]])
                    card += deg
                    max_deg = max(max_deg, deg)

        if card != cat2[v_list + ',' + label_seq]:
            print('ERROR: cardinality sanity check failed')
            print('   ' + str(cat2[v_list + ',' + label_seq]) + ' vs ' + str(card))
            return
        pcat_max_deg[sec[0] + '-' + sec[2] + ',' + sec[1] + ',' + first[0] + '-' + first[2] + ',' + first[1]] = max_deg


def persist(cat_file, pcat_file, pcat_max_deg_file):
    with open(cat_file, 'w') as f:
        for entry in cat2:
            f.write(template_id + ',' + entry + ',' + str(cat2[entry]) + '\n')
        for entry in cat1:
            f.write(template_id + ',' + entry + ',' + str(cat1[entry]) + '\n')

    prefix2 = ';'.join(['1-1'] * query_len) + ',0-0;0-0'
    prefix1 = ';'.join(['1-1'] * query_len) + ',0-0'
    with open(pcat_file, 'w') as f:
        for entry in cat2:
            f.write(template_id + ',' + prefix2 + ',' + entry + ',' + str(cat2[entry]) + '\n')
        for entry in cat1:
            f.write(template_id + ',' + prefix1 + ',' + entry + ',' + str(cat1[entry]) + '\n')

    with open(pcat_max_deg_file, 'w') as f:
        for entry in pcat_max_deg:
            f.write(template_id + ',' + prefix2 + ',' + entry + ',' + str(pcat_max_deg[entry]) + '\n')


def compute_cat1(query):
    for edge in query:
        entry = edge[0] + '-' + edge[2] + ',' + edge[1]
        if entry not in cat1:
            cat1[entry] = len(label2srcdest[edge[1]])

def read_graph(graph_file):
    with open(graph_file, 'r') as f:
        for line in f:
            e = line.strip().split(',')
            if e[1] not in label2srcdest:
                label2srcdest[e[1]] = []
            label2srcdest[e[1]].append((e[0], e[2]))

            if e[0] not in src2label2dst:
                src2label2dst[e[0]] = {}
            if e[1] not in src2label2dst[e[0]]:
                src2label2dst[e[0]][e[1]] = []
            src2label2dst[e[0]][e[1]].append(e[2])

            if e[2] not in dst2label2src:
                dst2label2src[e[2]] = {}
            if e[1] not in dst2label2src[e[2]]:
                dst2label2src[e[2]][e[1]] = []
            dst2label2src[e[2]][e[1]].append(e[0])

def read_queries(query_file):
    with open(query_file, 'r') as f:
        for line in f:
            e = line.strip().split(',')
            v_list = re.split(';|-', e[0])
            labels = e[1].split('->')
            query = []
            for i in range(0, len(v_list), 2):
                query.append((v_list[i], labels[int(i / 2)], v_list[i + 1]))
            queries.append(query)

def decompose(query):
    decoms = set([])
    for i in range(len(query)):
        for j in range(i + 1, len(query)):
            e1 = query[i]
            e2 = query[j]
            if e1[0] == e2[0] or e1[0] == e2[2] or e1[2] == e2[0] or e1[2] == e2[2]:
                decoms.add((e1, e2))
    return decoms

if __name__ == '__main__':
    main()
