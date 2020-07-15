import sys
import re


def main():
    output = set([])
    with open(sys.argv[1], 'r') as cat:
        for line in cat:
            entry = line.strip().split(',')
            v_list = re.split('-|;', entry[1])
            labels = entry[2].split('->')
            enc = encode(v_list, labels)
            output.add(enc + ',' + entry[3])

    with open(sys.argv[2], 'w') as out:
        for line in output:
            out.write(line + '\n')


def encode(v_list, labels):
    if len(labels) == 1:
        return '0,' + labels[0]
    elif len(labels) == 2:
        v0 = int(v_list[0])
        v1 = int(v_list[1])
        v2 = int(v_list[2])
        v3 = int(v_list[3])
        l0 = int(labels[0])
        l1 = int(labels[1])

        if v0 == v2 or v0 == v3:
            if v0 == v2:
                if l0 < l1:
                    return '0;0,' + labels[0] + ';' + labels[1]
                else:
                    return '0;0,' + labels[1] + ';' + labels[0]
            else:
                return '0;1,' + labels[0] + ';' + labels[1]
        else:
            if v1 == v2:
                return '0;1,' + labels[1] + ';' + labels[0]
            else:
                if l0 < l1:
                    return '1;1,' + labels[0] + ';' + labels[1]
                else:
                    return '1;1,' + labels[1] + ';' + labels[0]



if __name__ == '__main__':
    main()
