
def main():
    dosages = open('../data/ukbb_files/ukb_hla_v2.txt', 'r')
    hap_names = dosages.readline().split()
    dosages.close()

    new_names_file = open('output/haplotype_names/ukb_to_asterisk_names.txt', 'w')
    new_names_file.write('ukb_names	literature_names\n')

    for name in hap_names:
        halves = name.split('_')
        if len(halves[1]) == 3:
            halves[1] = '0' + halves[1]
        new_name = halves[0] + '*' + halves[1][:2] + ':' + halves[1][2:]
        new_names_file.write(name + '	' + new_name + '\n')
    new_names_file.close()
main()
