#!/usr/bin/env python3

import csv
import os
import shutil
import subprocess
from pathlib import Path

import fire
import xlsxwriter


def export_genotypes_xlsx_from_csvs(csv_dir, xlsx_output_path,
        csv_path_good='latest.good.reportTable.csv',
        csv_path_bad='latest.bad.reportTable.csv', 
        csv_path_genoset='latest.genoset.reportTable.csv',
    ):
    '''
    Aggregate 3 csv files into an xlsx file
        - 'latest.good.reportTable.csv',
        - 'latest.bad.reportTable.csv',
        - 'latest.genoset.reportTable.csv'

    Args:
        csv_dir: str
            xx
        xlsx_output_path: str
            xx
    Outputs:
        xlsx_file: file
            xx
    '''

    if csv_dir is None:
        csv_dir = ''

    csv_dir = Path(csv_dir)
    xlsx_output_path = Path(xlsx_output_path)
    xlsx_output_path.parent.mkdir(parents=True,exist_ok=True)

    workbook = xlsxwriter.Workbook(str(xlsx_output_path))

    workbook.set_properties({
        'title': 'Excel format Genome Report',
        'author': 'Jose Afonso Guerra-Assuncao',
        'comments': 'Generated for the Personal Genomes Project - United Kingdom Study',
    })

    input_csv_paths = {
        "Possibly Beneficial": csv_dir/csv_path_good,
        "Possibly Harmful": csv_dir/csv_path_bad,
        "Genosets": csv_dir/csv_path_genoset,
    }

    urlformat = workbook.add_format({'color': 'black', 'underline': 1})
    headerformat = workbook.add_format({'color': 'black', 'bold': 1})

    for ext in sorted(input_csv_paths.keys()):
        worksheet = workbook.add_worksheet(ext)
        cr = 0
        maxSize = []

        filename = input_csv_paths[ext]
        with open(filename, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                if not row:
                    continue
                for cc in range(len(row)):
                    data = row[cc]
                    strLength = 0
                    if 'href{' in data:
                        # url_start = data.index('href{') + 5
                        # url_end = data.index('}{')
                        # text_start = url_end + 2
                        # text_end = data.index('}')
                        # url = data[url_start:url_end]
                        # text = data[text_start:text_end]

                        url_parts = data.split("}{")
                        url = url_parts[0].replace("href{", "")
                        text = url_parts[1].replace("}", "")

                        worksheet.write_url(cr, cc, url, urlformat, text)
                        strLength = len(text)
                    else:
                        if cr == 0:
                            worksheet.write(cr, cc, data, headerformat)
                        else:
                            worksheet.write(cr, cc, data)
                        strLength = len(data)
                    if len(maxSize) > cc:
                        if strLength > maxSize[cc]:
                            maxSize[cc] = strLength
                    else:
                        maxSize.append(strLength)
                cr += 1

        for cc in range(len(maxSize)):
            worksheet.set_column(cc, cc, maxSize[cc])

    workbook.close()

def render_latex(tables_dir, sample_name, output_dir,
                 plot_path=None, vep_dir=None,
                 template_dir="templates", latex_template=None,
                 ):
    """
    Render the latex template with the tables and the version table

    Template
    -

    Latex render files:
        Common:
            latexTemplate.tex
            versionTable.txt
            GeneStructure.pdf

        Personal:
            SampleName.txt
            latest.good.reportTable.csv
            latest.bad.reportTable.csv
            latest.genoset.reportTable.csv

        with VEP:
            latest.summary.csv
            latest.polyphen_table.csv
            latest.var_class_table.csv
            latest.var_cons_table.csv

    Parameters
    ----------
    tables_dir
    sample_name
    output_dir

    Returns
    -------

    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    template_dir = Path(template_dir)

    if latex_template is None:
        if vep_dir is not None:
            latex_template = template_dir/f"reportTemplate_withVEP.tex"
        else:
            latex_template = template_dir/f"reportTemplate_ohneVEP.tex"

    # Prepare files for latex rendering
    ## Main file
    latex_main_name = f"{sample_name}_report"
    shutil.copy(latex_template, f"{output_dir}/{latex_main_name}.tex")

    ## Common files
    shutil.copy(template_dir/f"versionTable.txt", output_dir)
    shutil.copy(template_dir/f"GeneStructure.pdf", output_dir)

    ## Personal files
    with open(f"{output_dir}/SampleName.txt", "w") as f:
        f.write(sample_name)
    tables_dir = Path(tables_dir)
    shutil.copy(tables_dir/f"latest.good.reportTable.csv", output_dir)
    shutil.copy(tables_dir/f"latest.bad.reportTable.csv", output_dir)
    shutil.copy(tables_dir/f"latest.genoset.reportTable.csv", output_dir)
    if plot_path is not None:
        shutil.copy(plot_path, output_dir/'AncestryPlot.pdf')

    ## VEP files
    if vep_dir is not None:
        vep_dir = Path(vep_dir)
        shutil.copy(vep_dir/f"latest.summary.csv", output_dir)
        shutil.copy(vep_dir/f"latest.polyphen_table.csv", output_dir)
        shutil.copy(vep_dir/f"latest.var_class_table.csv", output_dir)
        shutil.copy(vep_dir/f"latest.var_cons_table.csv", output_dir)


    # Render latex
    cwd = os.getcwd()
    os.chdir(f"{output_dir}")
    for _ in range(3):
        cmd = f"pdflatex -interaction=nonstopmode {latex_main_name}.tex 2> /dev/null > /dev/null"
        ret = subprocess.run(cmd, shell=True)
    os.chdir(cwd)


funcs = {
    'export_genotypes_xlsx_from_csvs': export_genotypes_xlsx_from_csvs,
    'render_latex': render_latex,
}


class Tools(object):
    def __init__(self):
        super(Tools, self).__init__()


if __name__ == '__main__':
    for k, v in funcs.items():
        setattr(Tools, k, staticmethod(v))
    fire.Fire(Tools)


