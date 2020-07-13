import sys

lines = []

with open(sys.argv[1], 'r') as result:
    for line in result:
        if 'txt' in line:
            lines.append(line.strip())

lines.sort()

total_runtime = 0
num_queries = 0

for line in lines:
    info = line.split()
    total_runtime += float(info[2])
    num_queries += 1

print('Avg Time: ' + str(total_runtime / num_queries))
