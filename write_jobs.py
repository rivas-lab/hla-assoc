import subprocess
import sys
import time

def write_job(file_name, phe, out_name):
    job_file = open(file_name, 'w')
    job_file.write("#!/bin/bash\n#\n")
    job_file.write("#SBATCH --job-name=" + out_name + "\n")
    job_file.write("#SBATCH --output=job_output/" + out_name + ".%j.out\n")
    job_file.write("#SBATCH --error=job_output/" + out_name + ".%j.err\n")
    job_file.write("#SBATCH --time=5:00:00\n")
    job_file.write("#SBATCH --qos=normal\n")
    job_file.write("#SBATCH -p owners\n")
    job_file.write("#SBATCH --nodes=1\n")
    job_file.write("#SBATCH --mem=20Gb\n") 
    job_file.write("ml load R\n")
    job_file.write("Rscript reg_test.R /share/PI/mrivas/data/ukbb/phefiles/highconfidenceqc/" + phe + " " + out_name)
    job_file.close()
  
def submit_job(file_name):
    try:
        subprocess.check_call('sbatch {}'.format(file_name), shell=True)
        sys.stderr.write('{} submitted on sherlock2\n'.format(file_name))
        success = True
    except:
        sys.stderr.write('{} not submitted on sherlock2\n'.format(file_name))
        success = False

def choose_phes(num = 'all'):
    if num == 'all':
        map_file = open("/share/PI/mrivas/data/ukbb/phefiles/highconfidenceqc/highconfidenceqc_map.txt", "r")
        ids = []
        for line in map_file.readlines():
            ids.append(line.split()[0])
        return ids

    else:
        name_to_id = {}
        name_to_id_file = open("/share/PI/mrivas/data/ukbb/phefiles/highconfidenceqc/highconfidenceqc_map.txt", "r")
        for line in name_to_id_file.readlines():
            line = line.split(" ", 1)
            name_to_id[line[1][:-1]] = line[0]
        name_to_id_file.close()

        freq_to_id = {}

        freq_to_id_file = open("/share/PI/mrivas/data/ukbb/phefiles/highconfidenceqc/17_07_04_highconfidence_summary.csv", "r")
        info = freq_to_id_file.readlines()
        info = info[1:]
        for line in info:
            line = line.split(",")
            try:
                freq_to_id[int(line[1][:-1])] = name_to_id[line[0]]
            except:
                print(line[1][:-1])
                print(line[0])

        freq_to_id_file.close()
        keylist = list(freq_to_id.keys())
        keylist.sort()
        ids = []
#        print(keylist[-num:])
        for key in keylist[-num:]:
            ids.append(freq_to_id[key])
        return ids

def main():
#    write_job("run_reg.sh", "HC303.phe", "out")
#    submit_job("run_hello.sh")
    phe_ids = choose_phes()
    print(phe_ids)
    for phe in phe_ids:
        write_job("run_reg.sh", "HC" + phe + ".phe", "HC" + phe)
        submit_job("run_reg.sh")
        time.sleep(1)
main()
