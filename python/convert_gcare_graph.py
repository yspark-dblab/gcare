import sys

edges = set([])

with open(sys.argv[1], 'r') as gcare_data_file:
    for line in gcare_data_file:
        if line.startswith('e'):
            if line.strip() in edges:
                print(line.strip())
            edges.add(line.strip())

with open(sys.argv[2], 'w') as converted_file:
    for edge in edges:
        e = edge.split()
        src = e[1]
        dst = e[2]
        label = e[3]
        converted_file.write(src + ',' + label + ',' + dst + '\n')
