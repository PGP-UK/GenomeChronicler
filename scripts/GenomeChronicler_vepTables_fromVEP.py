#!/usr/bin/env python3

import subprocess
import sys
import re
from pathlib import Path

import fire

def get_vepTables_from_VEP(filename, out_dir="."):

    # # Sort out input
    # if len(sys.argv) != 3:
    #     sys.exit("Usage: {} VEPhtml OutputDir".format(sys.argv[0]))

    # filename = sys.argv[1]
    # out_dir = sys.argv[2]

    # if not os.path.exists(filename):
    #     sys.exit("Major error here (File not found)")

    out_file = open("{}/latest.summary.csv".format(out_dir), "w")
    out_file.write("Feature\tCount\n")

    recording = 0
    # with (filename, "r") as in_file:
    with subprocess.Popen(f'cat {filename} | grep "gen_stats"', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as proc:
        for line in proc.stdout:
            line = line.decode()
            # line = line.strip()
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

    with open(filename, "r") as in_file:
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

                out_file = open("{}/latest.{}.csv".format(out_dir, tab_name), "w")

                head = table.pop(0)
                head[0] = "Label"
                out_file.write("{}\n".format(",".join(head)))

                other = 0
                for d in table:
                    perc = f"{round(float(d[1]) / summer * 100, 2):0.2f}"
                    if float(perc) > 1:
                        d[0] = d[0].capitalize()
                        d[1] = perc
                        out_file.write("{}\n".format(",".join([str(x) for x in d])))
                    else:
                        other += int(d[1])

                if other > 0:
                    perc = f"{round(float(other) / summer * 100, 2):0.2f}"
                    out_file.write(f"Others,{perc}\n")
                    # out_file.write("{}\n".format(",".join([str(x) for x in d])))
                   


if __name__ == '__main__':
    fire.Fire(get_vepTables_from_VEP)
