import subprocess
import sys
import time

def write_processing_job(file_name, phe):
    job_name = phe + "proc"
    job_file = open(file_name, 'w')
    job_file.write("#!/bin/bash\n#\n")
    job_file.write("#SBATCH --job-name=" + job_name + "\n")
    job_file.write("#SBATCH --output=job_output/" + job_name + ".%j.out\n")
    job_file.write("#SBATCH --error=job_output/" + job_name + ".%j.err\n")
    job_file.write("#SBATCH --time=10:00\n")
    job_file.write("#SBATCH --qos=normal\n")
    job_file.write("#SBATCH -p owners\n")
    job_file.write("#SBATCH --nodes=1\n")
    job_file.write("#SBATCH --mem=2Gb\n") 
    job_file.write("source activate h2-estimation\n")
    job_file.write("Rscript examine_output.R reg_results/" + phe + "\n")
    job_file.write("python process_output.py " + phe)
    job_file.close()

def write_job_R(file_name, phe, out_name, phe_directory):
    job_name = out_name + "R"
    job_file = open(file_name, 'w')
    job_file.write("#!/bin/bash\n#\n")
    job_file.write("#SBATCH --job-name=" + job_name + "\n")
    job_file.write("#SBATCH --output=job_output/" + job_name + ".%j.out\n")
    job_file.write("#SBATCH --error=job_output/" + job_name + ".%j.err\n")
    job_file.write("#SBATCH --time=12:00:00\n")
    job_file.write("#SBATCH --qos=normal\n")
    job_file.write("#SBATCH -p owners\n")
    job_file.write("#SBATCH --nodes=1\n")
    job_file.write("#SBATCH --mem=20Gb\n") 
    job_file.write("ml load R\n")
    job_file.write("Rscript reg_test.R " + phe_directory + phe + " " + out_name)
    job_file.close()

def write_job_PLINK(file_name, phe, out_name, phe_directory):
    job_name = out_name + "_PLINK"
    ukbb_files = "$SCRATCH/ukbb_files/"
    job_file = open(file_name, 'w')
    job_file.write("#!/bin/bash\n#\n")
    job_file.write("#SBATCH --job-name=" + job_name + "\n")
    job_file.write("#SBATCH --output=job_output/" + job_name + ".%j.out\n")
    job_file.write("#SBATCH --error=job_output/" + job_name + ".%j.err\n")
    job_file.write("#SBATCH --time=1:00:00\n")
    job_file.write("#SBATCH --qos=normal\n")
    job_file.write("#SBATCH -p owners\n")
    job_file.write("#SBATCH --nodes=1\n")
    job_file.write("#SBATCH --mem=20Gb\n") 
#    job_file.write("ml load R\n")
#    job_file.write("source activate h2-estimation\n")
    job_file.write("ml load plink2\n")

    job_file.write("sampleqc=" + ukbb_files + "sampleqc_16698.phe\n")
    job_file.write("pheno=" + phe_directory + phe + "\n")

    job_file.write("bed=" + ukbb_files + "ukb_hla_v3.bed\nbim=" + ukbb_files + "ukb_hla_v3.bim\nfam=" + ukbb_files + "ukb_hla_v3.fam\n")
    job_file.write("out=$SCRATCH/PLINK_results/" + out_name + "\n")
    job_file.write("covar=" + ukbb_files + "covar_16698.phe\n")
    job_file.write("plink2 --memory 20000 --threads 4 --remove $sampleqc --bed $bed --bim $bim --fam $fam --out $out --pheno $pheno --covar $covar --covar-name age sex PC1 PC2 PC3 PC4 --glm firth-fallback")

    job_file.close()

  
def submit_job(file_name):
    try:
        subprocess.check_call('sbatch {}'.format(file_name), shell=True)
        sys.stderr.write('{} submitted on sherlock2\n'.format(file_name))
        success = True
    except:
        sys.stderr.write('{} not submitted on sherlock2\n'.format(file_name))
        success = False

def choose_phes_2(num = "all"):
    """Just list names, no intelligent sorting"""
    phes = []
    for i in range(446):
        phes.append(str(i))
    if num == "all":
        return phes
    else:
        return phes[:num]

def choose_phes(phe_directory, num = 'all'):
    if num == 'all':
        map_file = open(phe_directory + "highconfidenceqc_map.txt", "r")
        ids = []
        for line in map_file.readlines():
            ids.append(line.split()[0])
        return ids

    else:
        name_to_id = {}
        name_to_id_file = open(phe_directory + "highconfidenceqc_map.txt", "r")
        for line in name_to_id_file.readlines():
            line = line.split(" ", 1)
            name_to_id[line[1][:-1]] = line[0]
        name_to_id_file.close()

        freq_to_id = {}

        freq_to_id_file = open(phe_directory + "17_07_04_highconfidence_summary.csv", "r")
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

    num = "all"

    # choose which regressions to run
    if len(sys.argv) < 2:
        raise ValueError("Not enough arguments")
    elif sys.argv[1] not in ["PLINK", "R", "proc"]:
        raise ValueError("invalid argument") 
    else:
        reg = sys.argv[1]

    if len(sys.argv) == 3:
        num = int(sys.argv[2])

#    phe_directory = "/oak/stanford/groups/mrivas/ukbb/16698/phe/highconfidenceqc/"
    phe_directory = "/share/PI/mrivas/data/ukbb/phefiles/highconfidenceqc/"

    phe_ids = choose_phes(phe_directory, num)
#    phe_ids = choose_phes_2(num)

#    phe_ids = [phe for phe in phe_ids if phe > 275]

    #phe_ids = ['29', '65']
    print(phe_ids)
    if reg == "R":
        file_name = "scripts/run_reg_R.sh"
    elif reg == "PLINK":
        file_name = "scripts/run_reg_PLINK.sh"
#    if reg in ["PLINK", "R"]:
#        file_name = "run_reg.sh"
    elif reg == "proc":
        file_name = "scripts/process_output.sh"

    # run jobs
    for phe in phe_ids:
        if reg == "R":
            write_job_R(file_name, "HC" + phe + ".phe", "HC" + phe, phe_directory)
        elif reg == "PLINK":
            write_job_PLINK(file_name, "HC" + phe + ".phe", "HC" + phe, phe_directory)
        elif reg == "proc":
            write_processing_job(file_name, phe)

        submit_job(file_name)
        time.sleep(1)
main()
