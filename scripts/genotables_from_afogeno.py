import os
import subprocess
from pathlib import Path
import re

import sqlite3
import fire
import pandas as pd
from tqdm import tqdm

if __name__ == '__main__':
    import sys
    sys.path.append(os.path.abspath(os.getcwd()))
from scripts.utils import ucfirst


def processGenoset(genoset, genosets, genotypes, positiveGenosets):
    # example: logic = ((( $v[0] + $v[1] + $v[2] + $v[3] + $v[4] ) >= 1) and (!$v[5] and !$v[6] and !$v[7] and !$v[8] and !$v[9]))
    if genoset not in genosets:
        result = 0
        return result
    logic = genosets[genoset][0]
    logic = logic.replace('$', '')  # new for Python
    logic = logic.replace('!', 'not ')  # new for Python
    vars = genosets[genoset][1].split(",")
    v = {}
    g = {}

    # BIG FAT WARNING: DO SOMETHING CLEVER ABOUT THE iIDs OTHERWISE WE WILL BE IN TROUBLE

    # Fill in vars
    for var in vars:
        varComponents = var.split("=")
        # eval(varComponents[0] + " = evaluateGenotype(varComponents[1])")
        exec(varComponents[0][
             1:] + " = evaluateGenotype(varComponents[1], genoset, genosets, genotypes, positiveGenosets)")  # remove '$'

    # Evaluate Genoset
    # result = 0
    if logic == '':
        result = 0
    else:
        result = int(eval(logic))
    # exec("result = 1 if " + logic)

    # Put data back into genotypes for future use
    genotypes[genoset] = result

    # If positive, add to positiveGenosets
    if (result == 1):
        positiveGenosets[genoset] = 1
    return result


def evaluateGenotype(query, genoset, genosets, genotypes, positiveGenosets):
    result = 0

    if re.match(r"(.+)\((.+)\)", query):
        snp, geno = re.findall(r"(.+)\((.+)\)", query)[0]
        if snp not in genotypes:
            print(f"Some big problem here hum? [{snp}]")
            return 0
        if geno in genotypes[snp]:
            return 1
        else:
            return 0
    elif re.match(r"(gs.+)", query, re.IGNORECASE):
        genoset = re.findall(r"(gs.+)", query, re.IGNORECASE)[0]
        if genoset in genotypes:
            return genotypes[genoset]
        else:
            print(f"WARNING: Genoset {genoset} should have been previously computed...")
            return processGenoset(genoset, genosets, genotypes, positiveGenosets)
    else:
        raise Exception(f"I'm sorry Dave, I'm afraid I can't do that [{query}]")

    return result


def cleanSummaryString(row):
    # Removing commas from summary so we can use a CSV format for the LaTeX table later on.
    row = row.replace(",", ":")
    row = row.replace("&", "and")

    # Tidying up the summaries to appease my OCD
    # row = row.capitalize()
    row = ucfirst(row)

    # Sorry, but I'll have to cut on the info here...
    row = row[:47]
    if len(row) == 47:
        row = row + "..."

    return row


def generate_GnomAD_url(rsid, sth3, sth3_sql):
    rv3 = sth3.execute(sth3_sql, (rsid,))
    counter = 0
    for row in sth3.fetchall():
        counter += 1

    if counter:
        return "\\href{{http://gnomad.broadinstitute.org/awesome?query={0}}}{{Link}}".format(rsid)
    else:
        return ""


def generate_get_evidence_url(rsid, sth4, sth4_sql):
    rv4 = sth4.execute(sth4_sql, (rsid,))
    counter = 0
    for row in sth4.fetchall():
        counter += 1

    if counter:
        return "\\href{{http://evidence.pgp-hms.org/{0}}}{{Link}}".format(rsid)
    else:
        return ""


def generate_SNPedia_url(rsid):
    # return "\\href{{https://www.snpedia.com/index.php/{0}}}{{{1}}}".format(rsid.capitalize(), rsid)
    return "\\href{{https://www.snpedia.com/index.php/{0}}}{{{1}}}".format(ucfirst(rsid), rsid)


def generate_Summary_url(rsid, summary):
    # return "\\href{{https://www.snpedia.com/index.php/{0}}}{{{1}}}".format(rsid.capitalize(), summary)
    return "\\href{{https://www.snpedia.com/index.php/{0}}}{{{1}}}".format(ucfirst(rsid), summary)


def generate_ClinVar_url(rsid, sth5, sth5_sql):
    rv5 = sth5.execute(sth5_sql, (rsid,))
    counter = 0
    for row in sth5.fetchall():
        counter += 1

    if counter:
        return "\\href{{http://www.ncbi.nlm.nih.gov/clinvar/?term={0}}}{{Link}}".format(rsid)
    else:
        return ""


def geno_tables_from_afogeno(afogeno_path, output_dir, verbose=False, version="22-146", reference_dir="reference"):
    """
    This function generates the geno tables from the afogeno output.

    Inputs:
        afogeno_path: path to the afogeno output
        reference_dir: path to the reference directory
            snpedia.db
            gnomad.db
            getevidence.db
            clinvar.db
            genosetDependencies.txt
            parsedGenosets.txt
    Outputs:
        latest.good.reportTable.csv"
        latest.bad.reportTable.csv
        latest.genoset.reportTable.csv

    Parameters
    ----------
    afogeno_path
    output_dir
    verbose
    version
    reference_dir

    Returns
    -------

    """
    # dir = ""
    # if "SINGULARITY_NAME" in os.environ:
    #     dir = "/GenomeChronicler/"

    reference_dir = Path(reference_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = output_dir/"temp"
    temp_dir.mkdir(parents=True, exist_ok=True)

    # Note: It is quite clear this code can be optimised
    tmpBlacklist = {"rs10156191": 1, "rs12094543": 1, "rs1667255": 1, "rs17672135": 1, "rs6445975": 1, "rs7659604": 1}

    print(f"This is script [{__file__}] part of the report generator pipeline version [{version}];\n")

    # Connect to databases
    snpedia_db = sqlite3.connect(reference_dir/f"snpedia.db")
    print("Opened database successfully [SNPedia]")

    gnomad_db = sqlite3.connect(reference_dir/f"gnomad.db")
    print("Opened database successfully [GnomAD]")

    getevidence_db = sqlite3.connect(reference_dir/f"getevidence.db")
    print("Opened database successfully [GetEvidence]")

    clinvar_db = sqlite3.connect(reference_dir/f"clinvar.db")
    print("Opened database successfully [ClinVar]")

    # Prepare statements
    sth1c = snpedia_db.cursor()
    # sth1c.execute('SELECT * FROM data WHERE data.id=?')
    sth1c_sql = 'SELECT * FROM data WHERE data.id=?'
    sth1f = snpedia_db.cursor()
    # sth1f.execute('SELECT * FROM flagged WHERE id=?')
    sth1f_sql = 'SELECT * FROM flagged WHERE id=?'
    sth1s = snpedia_db.cursor()
    # sth1s.execute('SELECT * FROM strand WHERE id=?')
    sth1s_sql = 'SELECT * FROM strand WHERE id=?'
    sth1g = snpedia_db.cursor()
    # sth1g.execute('SELECT * FROM genoset WHERE id=?')
    sth1g_sql = 'SELECT * FROM genoset WHERE id=?'

    sth3 = gnomad_db.cursor()
    # sth3.execute('SELECT * FROM data WHERE rsid=?')
    sth3_sql = 'SELECT * FROM data WHERE rsid=?'

    sth4 = getevidence_db.cursor()
    # sth4.execute('SELECT * FROM data WHERE dbsnp_id=?')
    sth4_sql = 'SELECT * FROM data WHERE dbsnp_id=?'

    sth5 = clinvar_db.cursor()
    # sth5.execute('SELECT * FROM data WHERE rsid=?')
    sth5_sql = 'SELECT * FROM data WHERE rsid=?'

    # Input Genoset Data
    allGenosets = []
    with open(reference_dir/f"genosetDependencies.txt", "r") as genoset_deps_file:
        for line in genoset_deps_file:
            line = line.strip()
            if not line:
                continue
            allGenosets.append(line.split("\t")[0])

    genosets = {}
    genotypes = {}
    with open(reference_dir/f"parsedGenosets.txt", "r") as parsed_genosets_file:
        for line in parsed_genosets_file:
            line = line.strip()
            if not line:
                continue
            name, *data = line.split("\t")
            name = name.lower()
            genosets[name] = data

    # output_path = f'{outdir}/latest.good.reportTable.csv'
    # if not Path(output_path).exists():
    if True:
        IN = open(afogeno_path)
        rows_GOOD = []
        rows_BAD = []
        Path(f"{temp_dir}/latest.good.reportTable.csv").write_text("Mag.,Identifier,Genotype,Summary,GnomAD,GetEvidence,ClinVar\n")
        Path(f"{temp_dir}/latest.bad.reportTable.csv").write_text("Mag.,Identifier,Genotype,Summary,GnomAD,GetEvidence,ClinVar\n")
        proc_good = subprocess.Popen(f"sort | uniq | sort -t \',\' -k1,1nr >> {temp_dir}/latest.good.reportTable.csv",
                                     shell=True, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        proc_bad = subprocess.Popen(f"sort | uniq | sort -t \',\' -k1,1nr >> {temp_dir}/latest.bad.reportTable.csv",
                                    shell=True, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    
        for line in tqdm(IN, disable=not verbose):
            # for line in IN:
            line = line.strip()
            if line == "":
                continue
            if line.startswith("#"):
                continue

            debugBuffer = ""

            data = line.split("\t")
            if len(data) < 12:
                data.append("")

            debugBuffer += "\n+++\n" + "\t".join(data) + "\n"

            strand = "plus"
            snp = data[3]

            # rvS = sth1s.execute(snp)
            rvS = sth1s.execute(sth1s_sql, (snp,))
            counter = 0
            for row in sth1s.fetchall():
                strand = row[1]
                counter += 1
            if counter > 1:
                raise Exception("Oh my Gawd!!! DNA has two strands (and it isn't a good thing in this case)")

            # rvFlag = sth1f.execute(snp)
            rvS = sth1f.execute(sth1f_sql, (snp,))
            flagID = 0
            for row in sth1f.fetchall():
                flagID = 1
            if flagID == 1 and (data[8] == "0/0" or data[8] == "./."):
                continue

            # strand
            if strand == "minus":
                data[9] = data[9].translate(str.maketrans("ACTGactg", "TGACtgac"))
                data[10] = data[10].translate(str.maketrans("ACTGactg", "TGACtgac"))

            genotype = [data[9], data[10]]
            extra = data[11]

            genotypes[snp] = {}
            if extra != "":
                genotypes[snp][extra] = 1
            genotypes[snp][genotype[0]] = 1
            genotypes[snp][genotype[1]] = 1
            genotypes[snp][";".join(genotype)] = 1

            try:
                rv = sth1c.execute(sth1c_sql, (snp,))
            except Exception as e:
                print(data)
                print(e)
                continue
            alleles = {}
            countGenotype = 0
            for row in sth1c.fetchall():
                debugBuffer += "ID = " + str(row[0]) + "\n"
                debugBuffer += "REPUTE = " + str(row[1]) + "\n"
                debugBuffer += "MAGNITUDE = " + str(row[2]) + "\n"
                debugBuffer += "ALLELE1 =  " + str(row[3]) + "\n"
                debugBuffer += "ALLELE2 =  " + str(row[4]) + "\n"
                debugBuffer += "SUMMARY =  " + str(row[5]) + "\n"
                debugBuffer += "RSID =  " + str(row[6]) + "\n"
                debugBuffer += "IID =  " + str(row[7]) + "\n"

                row = list(row)
                row[5] = cleanSummaryString(row[5])

                a1 = row[3]
                a2 = row[4]

                # if row[2] is None or row[2] == 0:
                if row[2] is None or row[2] == 0 or row[2] == "0": # Excluding magnitude 0
                    continue

                alleles.setdefault(a1, {})[a2] = [row[1], row[0], row[2], f"({a1};{a2})", row[5]]
                alleles.setdefault(a2, {})[a1] = alleles[a1][a2]

                countGenotype += 1



            if countGenotype:
                if genotype[0] in alleles and genotype[1] in alleles[genotype[0]]:

                    geno_desc = alleles[genotype[0]][genotype[1]][4]
                    #Implement some filters here
                    if re.findall(r"^\s*$|^Common|^None|^Normal|^Average|^Benign most likely|^Ancestral value|^Unaffected Genotype", geno_desc, re.IGNORECASE):
                        continue
                    if re.findall(r"^Benign variant|^L22- S142-|^Depends on|^Extensive metabolizer", geno_desc, re.IGNORECASE):
                        continue
                    if re.findall(r"^Typical BuChE|Complex; generally normal risk|common in clinvar|Major allele, normal risk|Slight if any|Likely to be benign|Most common genotype", geno_desc, re.IGNORECASE):
                        continue
                    if re.findall(r"^Most likely a benign polymorphism|Benign \(harmless\) variant|^1.3x risk$", geno_desc, re.IGNORECASE):
                        continue
                    if re.findall(r"normal risk of migraine|More likely to go bald; common|Most likely benign polymorphism|Slightly increased lifespan?|^Benign polymorphism", geno_desc, re.IGNORECASE):
                        continue
                    if re.findall(r"No increased risk of |Likely to be a benign variant|Carrier of a benign change|Classified as benign in ClinVar", geno_desc, re.IGNORECASE):
                        continue

                    rsid = alleles[genotype[0]][genotype[1]][1]
                    gnomAD = generate_GnomAD_url(rsid, sth3, sth3_sql)
                    getevidence = generate_get_evidence_url(rsid, sth4, sth4_sql)
                    snpedia = generate_SNPedia_url(rsid)
                    clinvar = generate_ClinVar_url(rsid, sth5, sth5_sql)

                    if gnomAD == "" and getevidence == "" and clinvar == "":
                        continue

                    if rsid in tmpBlacklist:
                        continue

                    alleles[genotype[0]][genotype[1]][1] = snpedia
                    alleles[genotype[0]][genotype[1]][4] = generate_Summary_url(rsid,
                                                                                alleles[genotype[0]][genotype[1]][
                                                                                    4])
                    tmpString = ",".join(
                        [alleles[genotype[0]][genotype[1]][2], alleles[genotype[0]][genotype[1]][1],
                         alleles[genotype[0]][genotype[1]][3], alleles[genotype[0]][genotype[1]][4]])
                    tmpString = tmpString.encode('ascii', 'ignore').decode()  # Removing non-ASCII characters
                    if alleles[genotype[0]][genotype[1]][2] is None or alleles[genotype[0]][genotype[1]][2] == "":
                        continue
                    elif alleles[genotype[0]][genotype[1]][0] is None or alleles[genotype[0]][genotype[1]][0] == "":
                        # print("WARNING: Go here and with a fair sense of justice determine is this is good or bad\n", tmpString, "\n")
                        continue
                    elif alleles[genotype[0]][genotype[1]][0] == "Good":
                        proc_good.stdin.write(f"{tmpString},{gnomAD},{getevidence},{clinvar}\n".encode())
                    elif alleles[genotype[0]][genotype[1]][0] == "Bad":
                        proc_bad.stdin.write(f"{tmpString},{gnomAD},{getevidence},{clinvar}\n".encode())

                elif countGenotype > 2:
                    print("IMPORTANT: ", debugBuffer)
                    print("IMPORTANT: ", debugBuffer)
                    print("IMPORTANT: Please debug this here as I couldn't find the right alleles [", strand, "] [",
                          countGenotype, "]\n\n\n")

        proc_good.stdin.close()
        proc_bad.stdin.close()


    positiveGenosets = {}
    for genoset in tqdm(allGenosets, disable=not verbose):
        # result = processGenoset(genoset)
        result = processGenoset(genoset, genosets, genotypes, positiveGenosets)
        # if result:
        #     positiveGenosets[result[0]] = result[1]

    print("\t +++ POSITIVES:\n", positiveGenosets)

    with open(f"{output_dir}/latest.genoset.reportTable.csv", "w") as geno:
        geno.write("Magnitude,Identifier,Summary\n")

    sort_process = subprocess.Popen(f"sort -t \",\" -k 1,1nr >> {output_dir}/latest.genoset.reportTable.csv", stdin=subprocess.PIPE, shell=True)
    # proc_bad = subprocess.Popen(f"sort | uniq | sort -t \',\' -k1,1nr >> {outdir}/latest.bad.reportTable.csv", 
    #             shell=True, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    for geno in sorted(positiveGenosets.keys()):
        rvG = sth1g.execute(sth1g_sql, (geno,))
        while True:
            row = sth1g.fetchone()
            if row is None:
                break
            if row[2] != "0":
                sort_process.stdin.write(f"{row[2]},{generate_SNPedia_url(row[0])},{generate_Summary_url(row[0], cleanSummaryString(row[3]))}\n".encode())

    sort_process.stdin.close()
    # sort_process.wait()

    print("Operation done successfully")
    snpedia_db.close()

    filter_final_report_csv(f"{temp_dir}/latest.good.reportTable.csv", f"{output_dir}/latest.good.reportTable.csv")
    filter_final_report_csv(f"{temp_dir}/latest.bad.reportTable.csv", f"{output_dir}/latest.bad.reportTable.csv")


def filter_final_report_csv(csv_input_path, csv_output_path=None, DEBUG=False):
    '''

    '''
    # read csv
    df = pd.read_csv(csv_input_path,dtype={'Mag.':'str'})

    # filter
    if len(df.columns) <= 4:
        if DEBUG:
            print("WARNING: No links found")
            return
    df = df[~(df['Mag.'].astype('float')==0)]
    df = df[~(df['Identifier'].str.startswith('{rs')|df['Identifier'].str.startswith('{i'))]

    # export
    if csv_output_path is None: # inplace
        # warning, inplace operation
        csv_output_path = csv_input_path
    Path(csv_output_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv_output_path, index=False)



funcs = {
    'geno_tables_from_afogeno': geno_tables_from_afogeno,
    'filter_final_report_csv': filter_final_report_csv,
}


class Tools(object):
    def __init__(self):
        super(Tools, self).__init__()


if __name__ == '__main__':
    for k, v in funcs.items():
        setattr(Tools, k, staticmethod(v))
    fire.Fire(Tools)
