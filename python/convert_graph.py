import sys

vertices = {} # map from old to new ID
current_vid_max = -1
edges = []
min_label = -1

with open(sys.argv[1], 'r') as data:
    for line in data:
        edge = line.strip().split(',')
        src = int(edge[0])
        label = int(edge[1])
        dest = int(edge[2])

        min_label = min(min_label, label)

        if src not in vertices:
            current_vid_max += 1
            vertices[src] = current_vid_max
        if dest not in vertices:
            current_vid_max += 1
            vertices[dest] = current_vid_max

        edges.append((str(vertices[src]), label, str(vertices[dest])))


sorted_vids = sorted(list(vertices.values()))
print('Min_EL: ' + str(min_label))

with open(sys.argv[2], 'w') as gcare_graph:
    gcare_graph.write('t # 0\n')
    for v in sorted_vids:
        gcare_graph.write('v ' + str(v) + ' 0\n')
    for e in edges:
        gcare_graph.write('e ' + e[0] + ' ' + e[2] + ' ' + str(e[1] - min_label) + '\n')

