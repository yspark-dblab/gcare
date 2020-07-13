import sys

lines = []

with open(sys.argv[1], 'r') as result:
    for line in result:
        if 'txt' in line:
            lines.append(line.strip())

lines.sort()

for line in lines:
    info = line.split()
    template_and_query_id = info[0].split('/')[-1]
    template_id = template_and_query_id.split('_')[0][1:]
    query_id = template_and_query_id.split('_')[1].split('.')[0]
    print(template_id + ',' + query_id + ',' + info[1])
