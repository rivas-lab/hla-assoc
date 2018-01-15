import pandas as pd
import numpy as np
#import rpy2.robjects as robjects
import sys

def main():
#    print("here")
#    phe = sys.argv[1]
    phe_list = []
#    cancer_phes = open("../data/cancerphes.txt","r")
#    for line in cancer_phes.readlines():
#        phe_list.append(line[:-1])
#    cancer_phes.close()
    
    for i in range(446):
        if i not in [29,65]:
            phe_list.append("HC" + str(i))

    for phe in phe_list:
        print("phe: " + phe)
        #plink_name = "hla_celiac.PHENO1.glm.logistic.hybrid"
    #    plink_name = "/scratch/users/jolivier/PLINK_results/" + phe + ".PHENO1.glm.logistic.hybrid"
        plink_name = "../data/PLINK_results/" + phe + ".PHENO1.glm.logistic.hybrid"
    
        R_name = "output/reg_test/reg_results/" + phe + "_rounded_add.txt"
        df = pd.read_csv(plink_name, delimiter='\t', header=0, usecols = [2, 5, 6, 8, 11])
    
        print("read in plink")
    
        # don't include covariate lines
        df = df.loc[df['TEST'] == 'ADD']
        df = df.set_index("ID")
    
        # initialize lines of dataframe to be filled in
        df["LOG OR (R)"] = "NA"
        df["OR (R)"] = "NA"
        df["P (R)"] = "NA" 
    
    #    R_name = "out_dosage_add.txt"
    
        # extract p values and odds ratios from R results
        R_file = open(R_name, 'r')
        hap = ''
        for line in R_file.readlines():
            if line[0] == '$':
                hap = line[1:-1]
    #            print(hap)
            else:
                line = line.split()
                if len(line) > 4:
                    if line[0] == "as.numeric(gcounts)":
                        df.at[hap, "LOG OR (R)"] = line[1]
                        if line[1] != 'NA':
                            df.at[hap, "OR (R)"] = np.exp(float(line[1]))
                        if line[4] == '<':
                            df.at[hap, "P (R)"] = line[4] + line[5]
                        else:
                            df.at[hap, "P (R)"] = line[4]
    
        print("read in R")
    
        df = df[df["OR (R)"] != "NA"]
        pd.to_numeric(df["OR (R)"])
    
        df["LOG OR"] = np.log(df["OR"])
    
        df["OR DIFF"] = np.power(df["OR (R)"] - df["OR"], 2)
     #   df = df.loc[df["FIRTH?"] == 'N']
        df = df.drop(["FIRTH?", "TEST"], axis=1)
    
        df.to_csv("output/process_output/" + phe + "_processed_firth.txt", sep="\t")
        print("saved csv")
#    df["DIFF > 0.001"] = df["OR DIFF"] > 0.001
##    df = df.loc[df["OR DIFF"] > 0.001 == True]
#    df = df.loc[df["DIFF > 0.001"] == True]
#    print(df)
#    plt.scatter(df["LOG OR"], df["LOG OR (R)"])
#    plt.savefig("out_dosage_add_plot.png")
main()
