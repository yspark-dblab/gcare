import sys
import os

def main():
    file_prefix = sys.argv[2]
    dir_name = '/'.join(file_prefix.split('/')[:-1])
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    queries = []
    with open(sys.argv[1], 'r') as data:
        for line in data:
            query = line.strip().split(',')
            v_list = query[0]
            label_seq = query[1]
            queries.append((v_list, label_seq))

    for i in range(len(queries)):
        create_query(queries[i], file_prefix, i + 1)


def create_query(query, dest_prefix, index):
    v_list_str = query[0]
    label_seq_str = query[1]
    v_list_edges = v_list_str.split(';')
    label_seq = label_seq_str.split('->')

    vertices = set([])
    edges = []
    for i in range(len(v_list_edges)):
        src_dest = v_list_edges[i].split('-')
        src = int(src_dest[0])
        dest = int(src_dest[1])
        vertices.add(src)
        vertices.add(dest)
        edges.append((str(src), label_seq[i], str(dest)))

    vertices = sorted(list(vertices))

    file_name = dest_prefix + '_' + '{:02d}'.format(index) + '.txt'
    with open(file_name, 'w') as query_file:
        query_file.write('t # s ' + str(index) + '\n')
        for v in vertices:
            query_file.write('v ' + str(v) + ' -1 -1\n')

        for e in edges:
            query_file.write('e ' + e[0] + ' ' + e[2] + ' ' + e[1] + '\n')



if __name__ == '__main__':
    main()
