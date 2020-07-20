import sys

queries = set([])
if len(sys.argv) == 4:
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    for q in range(start, end + 1):
        queries.add('Q' + str(q))

lines = []

with open(sys.argv[1], 'r') as result:
    for line in result:
        if 'txt' in line:
            if len(queries) == 0 or line.split()[0].split('/')[-2] in queries:
                lines.append(line.strip())

lines.sort()

total_runtime = 0
num_queries = 0

for line in lines:
    info = line.split()
    total_runtime += float(info[2])
    num_queries += 1

print('Avg Time: ' + str(total_runtime / num_queries))
