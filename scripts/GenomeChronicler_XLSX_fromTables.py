#!/usr/bin/env python3

# import sys
# import os
import csv
from pathlib import Path

# import pandas as pd
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

if __name__ == '__main__':
    fire.Fire(export_genotypes_xlsx_from_csvs)
