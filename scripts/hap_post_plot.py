#import matplotlib.pyplot as plt
import pandas as pd

def list_to_float(a):
    new_a = []
    for l in a:
       l = l.replace('"', "")
       l = l.replace("'", "")
       try:
           new_a.append(float(l))
       except:
          pass 
    return new_a

def dot(K, L):
    dot_res = 0
    for i in range(len(K)):
        dot_res += K[i] * L[i]
        #print("K: {} L: {} dot: {}".format(K[i], L[i], dot_res))
    return dot_res   


def line_quotes(string, num = True):
    """Get all the items within quotes in a line"""
    flip = False
    curr_word = ""
    words = []
    for i in range(len(string)):
        if string[i] == '"' and flip:
            words.append(curr_word)
            curr_word = ""
            flip = False
        elif string[i] == '"' and not flip:
            flip = True
        elif flip:
            curr_word += string[i]
    #print("words: {}".format(words))
    new_words = []
    for w in words:
        try:
            float(w)
            # if num, save number; otherwise append 1 
            if num:
                new_words.append(float(w))
            else:
                new_words.append(1.)
        # append 0 if it's not a string
        except:
            new_words.append(0) 
    return new_words

def parse_bma(fn, phe, df, df_EV, df_SD):
    """Parse bma file, getting posterior probabilities, postmean, and postse"""
    f = open(fn, "r")
    rec = False
    row = 0
    hap_probs = {}
    post = []
    for line in f.readlines():
        words = line.split()
        if words[0] == '""':
            rec = False
        if rec:
            hap = words[0].split(".")[0]
            if row > 1:
                hap_probs[hap] = hap_probs[hap] + line_quotes(line, num=False)
            else:
                sep_line = line_quotes(line, num=True)
                hap_probs[hap] = line_quotes(line, num=False)[3:]
                #print("column values: {}".format(df_EV.columns.values))
                #print("allele: {}".format(hap))
                #print(df_EV.columns.values)
		try: 
                    df_EV[hap][phe] = sep_line[1]
                    df_SD[hap][phe] = sep_line[2]
		except:
                    print("failed: {}, {}".format(phe, hap))
        if words[0] == "Intercept":
            rec = True
            row += 1
        if words[0] == "post":
#            print("post: ")
#            print(list_to_float(words))
      #      print("words: {}".format(words))
            post += list_to_float(words)
      #      print("post: {}".format(post))
    for hap in hap_probs.keys():
#        print(hap_probs[hap])
#        print("post: {}".format(post))
        df[hap][phe] = dot(hap_probs[hap], post)
#        print("entry: {}".format(df[hap][phe]))
#    print(post)
#    print(hap_probs)
    f.close()
    return df, df_EV, df_SD

def main():
    
    BMA_phes = open("../notebooks/output/compare_cutoff_pvals/BMA_phe_torun.txt","r")
    phes = BMA_phes.readline().split()
    print(len(phes))
#    phes.remove("HC322")
    print(len(phes))
    BMA_phes.close()

    for i in range(len(phes)):
        if phes[i][:2] != "HC":
            phes[i] = "cancer" + phes[i]
    #print("phes: {}".format(phes))
    print("number of phenotypes: {}".format(len(phes))) 

#    phes = []
#    # Create phe names for files we have 
#    for i in range(446):
#        if i not in [3,10, 20, 29, 65, 74, 107, 127, 168, 205, 232, 236, 291,317,340,372,405,410,426,427]:
#            phes.append("HC" + str(i))
#    cancerphes = open("../data/cancerphes.txt", "r")
#    for line in cancerphes.readlines():
#        if int(line[:-1]) not in [1006, 1087, 1045, 1042, 1040, 1064]:
##            phes.append(line[:-1])
#            phes.append("cancer" + line[:-1])
#
#            print(line[:-1])
#    cancerphes.close()

#    # read in haplotypes (with > 5 occurences)
#    hap_names = open("../data/curr_haps.txt", "r")
#    haps = hap_names.readline().split()
#    hap_names.close()

    BMA_alleles = open("../notebooks/output/compare_cutoff_pvals/BMA_allele_torun.txt","r")
    haps = BMA_alleles.readline().split()
    BMA_alleles.close()

#    print("alleles: {}".format(haps))
    print("number of allelotypes: {}".format(len(haps)))

    # initialize dataframes
    df = pd.DataFrame(index = phes, columns = haps, dtype=float)
    df_EV = pd.DataFrame(index = phes, columns = haps, dtype=float)
    df_SD = pd.DataFrame(index = phes, columns = haps, dtype=float)
    print("index: {}".format(df.index.values))

 #   print(df.index)

    # parse bma files and put relevant info into dfs
    for phe in phes:
        #print("phe: " + phe)
        if phe[:6] == "cancer":
            phe_name = phe[6:]
        else:
  	    phe_name = phe
        df, df_EV, df_SD = parse_bma("output/reg_test/results_BMA/bma_" + phe_name + "_adjp.txt", phe, df, df_EV, df_SD)
    df.to_csv("output/hap_post_plot/phe_hap_table_post_adjp.csv")
    df_EV.to_csv("output/hap_post_plot/phe_hap_table_EV_adjp.csv")
    df_SD.to_csv("output/hap_post_plot/phe_hap_table_SD_adjp.csv")

    genes = ['A', 'C', 'B', 'DQB1', 'DPB1', 'DPA1', 'DQA1', 'DRB4', 'DRB5', 'DRB1', 'DRB3']

    # create individual gene tables
    for gene in genes:
        filter_col = [col for col in df if col.startswith(gene)]
        df_gene = df[filter_col]
        df_gene.to_csv("output/hap_post_plot/phe_post_table_" + gene + "_adjp.csv")
#    for gene in genes:
#       filter_col = [col for col in df_EV if col.startswith(gene)]
#       df_gene = df[filter_col]
#       df_gene.to_csv("~/oak/users/jolivier/repos/hla-assoc/phe_post_table_" + gene + ".csv")
#    for gene in genes:
#       filter_col = [col for col in df if col.startswith(gene)]
#       df_gene = df[filter_col]
#       df_gene.to_csv("~/oak/users/jolivier/repos/hla-assoc/phe_post_table_" + gene + ".csv")
#
main()
