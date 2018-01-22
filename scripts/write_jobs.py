import subprocess
import sys
import time


def write_printrds_job(file_name, phe):
    job_name = phe + "print"
    job_file = open(file_name, 'w')
    job_file.write("#!/bin/bash\n#\n")
    job_file.write("#SBATCH --job-name=" + job_name + "\n")
    job_file.write("#SBATCH --output=output/write_jobs/" + job_name + ".%j.out\n")
    job_file.write("#SBATCH --error=output/write_jobs/" + job_name + ".%j.err\n")
    job_file.write("#SBATCH --time=10:00\n")
    job_file.write("#SBATCH --qos=normal\n")
    job_file.write("#SBATCH -p owners\n")
    job_file.write("#SBATCH --nodes=1\n")
    job_file.write("#SBATCH --mem=2Gb\n") 
    job_file.write("ml load R\n")
#    job_file.write("Rscript print_rds.R results_BMA/bma_" + phe + "_10.rds " + "results_BMA/bma_" + phe + "_10.txt")
    if phe in ["HC303","HC382","HC38","HC430","HC219"]:
        job_file.write("Rscript print_rds.R output/reg_test/results_BMA/bma_" + phe + "_round_all.rds " + "output/reg_test/results_BMA/bma_" + phe + "_round_adjp_all.txt")
    else:
        job_file.write("Rscript print_rds.R output/reg_test/results_BMA/bma_" + phe + "_10_round_all.rds " + "output/reg_test/results_BMA/bma_" + phe + "_round_adjp_all.txt")


    job_file.close()

def write_processing_job(file_name, phe):
    job_name = phe + "proc"
    job_file = open(file_name, 'w')
    job_file.write("#!/bin/bash\n#\n")
    job_file.write("#SBATCH --job-name=" + job_name + "\n")
    job_file.write("#SBATCH --output=output/write_jobs/" + job_name + ".%j.out\n")
    job_file.write("#SBATCH --error=output/write_jobs/" + job_name + ".%j.err\n")
    job_file.write("#SBATCH --time=10:00\n")
    job_file.write("#SBATCH --qos=normal\n")
    job_file.write("#SBATCH -p owners\n")
    job_file.write("#SBATCH --nodes=1\n")
    job_file.write("#SBATCH --mem=2Gb\n") 
    job_file.write("source activate h2-estimation\n")
    job_file.write("Rscript examine_output.R output/reg_test/reg_results/" + phe + "\n")
    job_file.write("python process_output.py " + phe)
    job_file.close()

# make sure to check the true/false settings to see what's being run
def write_job_R(file_name, phe, out_name, phe_directory):
    job_name = out_name + "R"
    job_file = open(file_name, 'w')
    job_file.write("#!/bin/bash\n#\n")
    job_file.write("#SBATCH --job-name=" + job_name + "\n")
    job_file.write("#SBATCH --output=output/write_jobs/" + job_name + ".%j.out\n")
    job_file.write("#SBATCH --error=output/write_jobs/" + job_name + ".%j.err\n")
    job_file.write("#SBATCH --time=24:00:00\n")
    job_file.write("#SBATCH --qos=normal\n")
    job_file.write("#SBATCH -p owners\n")
    job_file.write("#SBATCH --nodes=1\n")
    job_file.write("#SBATCH --mem=20Gb\n") 
#    job_file.write("ml load R\n")
    job_file.write("source activate h2-estimation\n")
    job_file.write("ml load libpng/1.2.57\n")
    job_file.write("ml load cairo/1.14.10\n")
    job_file.write("Rscript reg_test.R " + phe_directory + phe + " " + out_name + " 10")
    job_file.close()

def write_job_PLINK(file_name, phe, phe_directory):
    job_name = phe + "_PLINK"
    ukbb_files = "../data/ukbb_files/"
    job_file = open(file_name, 'w')
    job_file.write("#!/bin/bash\n#\n")
    job_file.write("#SBATCH --job-name=" + job_name + "\n")
    job_file.write("#SBATCH --output=output/write_jobs/" + job_name + ".%j.out\n")
    job_file.write("#SBATCH --error=output/write_jobs/" + job_name + ".%j.err\n")
    job_file.write("#SBATCH --time=1:00:00\n")
    job_file.write("#SBATCH --qos=normal\n")
    job_file.write("#SBATCH -p owners\n")
    job_file.write("#SBATCH --nodes=1\n")
    job_file.write("#SBATCH --mem=20Gb\n") 
#    job_file.write("ml load R\n")
#    job_file.write("source activate h2-estimation\n")
    job_file.write("ml load plink2\n")

    job_file.write("sampleqc=" + ukbb_files + "sampleqc_16698.phe\n")
    job_file.write("pheno=" + phe_directory + phe + ".phe\n")

    job_file.write("bed=" + ukbb_files + "ukb_hla_v3.bed\nbim=" + ukbb_files + "ukb_hla_v3.bim\nfam=output/create_fam/ukb_hla_v3.fam\n")
    job_file.write("out=data/PLINK_results/" + phe + "\n")
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

def choose_phes_HC_2(num = "all"):
    """Just list names, no intelligent sorting"""
    phes = []
    for i in range(446):
        phes.append("HC" + str(i))
    if num == "all":
        return phes
    else:
        return phes[:num]

def choose_phes_cancer(phe_directory, num = "all"):
    cancermap = open(phe_directory + "cancermap.txt", "r")
    phes = []
    for line in cancermap.readlines():
        phes.append(line.split()[0])
    cancermap.close()
    if num == "all":
        return phes
    else:
        return phes[:num]
        

def choose_phes_HC(phe_directory, num = 'all'):
    """return list of phe ids to use; no HC in front"""
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
    elif sys.argv[1] not in ["BMA","PLINK", "R", "proc", "print"]:
        raise ValueError("invalid argument") 
    else:
        reg = sys.argv[1]

    if len(sys.argv) == 3:
        num = int(sys.argv[2])

    # choose which phenotype set to run on
    phe_type = "HC"

#    phe_directory = "/oak/stanford/groups/mrivas/ukbb/16698/phe/highconfidenceqc/"
    if phe_type == "HC":
        phe_directory = "/share/PI/mrivas/data/ukbb/phefiles/highconfidenceqc/"
        phe_ids = choose_phes_HC_2(num)

    elif phe_type == "cancer":
       phe_directory = "/oak/stanford/groups/mrivas/private_data/ukbb/16698/phe/cancer3/"
       phe_ids = choose_phes_cancer(phe_directory, num)

    # only run on phenotypes we perform BME
    # analysis for
    BMA = False
#    phe_ids = [3,10, 20, 29, 65, 74, 107, 127, 168, 205, 232, 236, 291,317,340,372,405,410,426,427]
#    phe_ids = ["3","10"]
    if BMA:
        phe_file = open("../notebooks/output/compare_cutoff_pvals/sig_phe_" + phe_type + ".txt")
        if num == "all": 
            phe_ids = phe_file.readline().split()
        else:
            phe_ids = phe_file.readline().split()[:num]

    print(phe_ids)
    #phe_ids = ["HC303","HC382","HC38","HC430","HC219"]
    if reg == "R":
        file_name = "shell_scripts/run_reg_R.sh"
    elif reg == "PLINK":
        file_name = "shell_scripts/run_reg_PLINK.sh"
#    if reg in ["PLINK", "R"]:
#        file_name = "run_reg.sh"
    elif reg == "proc":
        file_name = "shell_scripts/process_output.sh"
    elif reg == "print":
        file_name = "shell_scripts/run_print.sh"
    elif reg == "BMA":
        file_name = "shell_scripts/run_BMA.sh"

    # run jobs
    for phe in phe_ids:
        if reg == "R" or reg == "BMA":
            write_job_R(file_name, phe + ".phe", phe, phe_directory)
        elif reg == "PLINK":
            write_job_PLINK(file_name, phe, phe_directory)
        elif reg == "proc":
            write_processing_job(file_name, phe)
        elif reg == "print":
            write_printrds_job(file_name, phe)

        submit_job(file_name)
        time.sleep(1)
main()
