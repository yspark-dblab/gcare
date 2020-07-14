import sys

def main():
    prev_template_id = ''
    lines = []
    
    with open(sys.argv[1], 'r') as converted_result:
        for line in converted_result:
            template_id = line.strip().split(',')[0]
            if template_id != prev_template_id and template_id != '201':
                if prev_template_id == '':
                    save(lines, '201')
                else:
                    save(lines, prev_template_id)
                lines = []
                prev_template_id = template_id
    
            lines.append(line.strip())

        save(lines, prev_template_id)

def save(lines, template_id):
    file_name = sys.argv[2] + '_' + template_id + '.csv'
    with open(file_name, 'w') as f:
        for line in lines:
            f.write(line + '\n')


if __name__ == '__main__':
    main()
        

        


