import re
import subprocess
from pathlib import Path

import fire

def ucfirst(s):
    if len(s) > 0:
        return s[0].upper() + s[1:]
    else:
        return s


def clean_bam_file_noCHR(bam_path, output_dir=None):
    """
    Reheader the BAM file and create an index
    - Remove the "chr" prefix from the chromosome names

    Input: .bam
    Output: .clean.bam, .clean.bam.bai

    Parameters
    ----------
    bam_path: str or Path
        BAM file path
    output_dir: str, Path or None
        - If None, the output file will be saved in the same directory as the input file
        - Else, the output file will be saved in the specified directory

    Returns
    -------
    None

    """
    # Reheader the BAM file and create an index
    bam_path = Path(bam_path)
    if output_dir is None:
        output_dir = bam_path.parent
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir/ f"{bam_path.stem}.clean.BAM"
    subprocess.run(
        f"samtools reheader -c 'perl -pe \"s/^(@SQ.*)(\tSN:)chr/\$1\$2/\"' {bam_path} > {output_path}",
        shell=True)

    print(f"\t +++ INFO: Indexing the BAM file")
    subprocess.run(f"samtools index -@ 6 {output_path}", shell=True)


def get_vep_tables_from_vep(vep_path, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    out_file = open(f"{output_dir}/latest.summary.csv", "w")
    out_file.write("Feature\tCount\n")

    recording = 0
    with subprocess.Popen(f'cat {vep_path} | grep "gen_stats"', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as proc:
        for line in proc.stdout:
            line = line.decode()
            line = line.rstrip("\n")
            if not line:
                continue

            # Note for future generations: This is as silly as Perl gets, but it does the job (for now). I will re-write this with proper parsing at some point, or even better, use an off-the-shelve HTML parsing module.
            line = line.replace("<tr>", "\n")
            line = re.sub(r"</table>", "\n</table>\n", line)
            line = re.sub(r"<td>", "\t", line)
            line = re.sub(r"</td>|</tr>", "", line)

            lines = line.split("\n")
            for l in lines:
                if "stats_table" in l:
                    recording += 1
                    continue
                elif ("</table>" in l.lower()) and (recording == 2):
                    recording = 0

                # l = l.strip()
                l = re.sub(r"^\s+", "", l)
                l = re.sub(r"\s+$", "", l)
                l = re.sub(r"\s{2,}", "\t", l)

                if ("variants processed" in l.lower()):
                    continue
                if ("lines of output written" in l.lower()):
                    continue



                if recording == 2:
                    out_file.write("{}\n".format(l))

    # in_file.close()
    out_file.close()

    with open(vep_path, "r") as in_file:
        for line in in_file:
            line = line.strip()
            if not line or "drawTable" not in line or "[" not in line:
                continue

            match = re.search(r"drawTable\('(.+?)'.+?\[(.+)\]", line)
            if match:
                tab_name = match.group(1)
                tab_string = match.group(2)

                counter = 0
                summer = 0
                table = []
                for pair in re.findall(r"\[(.+?)\]", tab_string):
                    pair = pair.replace("'", "")
                    pair = pair.replace("_", " ")
                    if counter == 0:
                        pair = pair.replace(" ", "")
                    split = pair.split(",")
                    table.append(split)

                    if counter > 0:
                        summer += int(split[1])

                    counter += 1

                out_file = open("{}/latest.{}.csv".format(output_dir, tab_name), "w")

                head = table.pop(0)
                head[0] = "Label"
                out_file.write("{}\n".format(",".join(head)))

                other = 0
                for d in table:
                    perc = f"{round(float(d[1]) / summer * 100, 2):0.2f}"
                    if float(perc) > 1:
                        d[0] = ucfirst(d[0])
                        d[1] = perc
                        out_file.write("{}\n".format(",".join([str(x) for x in d])))
                    else:
                        other += int(d[1])

                if other > 0:
                    perc = f"{round(float(other) / summer * 100, 2):0.2f}"
                    out_file.write(f"Others,{perc}\n")
                    # out_file.write("{}\n".format(",".join([str(x) for x in d])))


funcs = {
    'clean_bam_file_noCHR': clean_bam_file_noCHR,
    'get_vep_tables_from_vep': get_vep_tables_from_vep,
}


class UtilsTools(object):
    def __init__(self):
        super(UtilsTools, self).__init__()


if __name__ == '__main__':
    for k, v in funcs.items():
        setattr(UtilsTools, k, staticmethod(v))
    fire.Fire(UtilsTools)

