import sys

vertices = set([])
edges = []

with open(sys.argv[1], 'r') as data:
    for line in data:
        edge = line.strip().split(',')
        src = edge[0]
        label = edge[1]
        dest = edge[2]

        vertices.add(int(src))
        vertices.add(int(dest))
        edges.append((src, label, dest))


vertices = sorted(list(vertices))

with open(sys.argv[2], 'w') as gcare_graph:
    gcare_graph.write('t # 0\n')
    for v in vertices:
        gcare_graph.write('v ' + str(v) + ' 0\n')
    for e in edges:
        gcare_graph.write('e ' + e[0] + ' ' + e[2] + ' ' + e[1] + '\n')

