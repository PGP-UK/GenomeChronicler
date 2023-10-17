import json
import os
import subprocess
from pathlib import Path
import re

import sqlite3
import fire
import pandas as pd
from tqdm import tqdm


def link_warper_split(x):
    tokens = x.split('}{')
    link = tokens[0].split('\\href{')[-1]
    text = tokens[1].strip('{}')
    return text, link

def row2entry(row):
    data = {}
    if 'Magnitude' in row:
        magnitude = row['Magnitude']
    elif 'Mag.' in row:
        magnitude = row['Mag.']
    else:
        magnitude = None
    data['magnitude'] = magnitude

    if 'Genotype' in row:
        data['genotype'] = row['Genotype'].strip('()').split(';')

    text, _ = link_warper_split(row['Summary'])
    desc_full = text
    desc_short = text
    if len(text) > 50:
        desc_short = text[:47] + '...'
    data['description'] = {
        'full': desc_full,
        'short': desc_short,
    }
    id_name, id_link = link_warper_split(row['Identifier'])
    data['identifier'] = {
        'link': id_link,
        'name': id_name,
    }
    db_links = []
    db_links += [{'link': id_link, 'name': 'SNPedia'}]
    for db_name in 'GnomAD GetEvidence ClinVar'.split():
        if db_name not in row:
            continue
        try:
            _, link = link_warper_split(row[db_name])
            db_links += [{'link': link, 'name': db_name}]
        except:
            db_links += [{'link': '', 'name': db_name}]
    data['db_links'] = db_links
    return data



def geno_tables_to_jsons(table_full_dir, output_dir):
    table_full_dir = Path(table_full_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    tables = {
        'beneficial_data.json': 'latest.good.reportTable.csv',
        'harmful_data.json': 'latest.bad.reportTable.csv',
        'genosets_data.json': 'latest.genoset.reportTable.csv',
    }
    for json_name, table_name in tables.items():
        path_table = table_full_dir / table_name
        df = pd.read_csv(path_table)
        data = df.apply(row2entry,axis=1).tolist()
        # Path()
        # data = table_to_json_data(df)
        path_json = output_dir / json_name
        path_json.write_text(json.dumps(data, indent=2))
        # pd.DataFrame(data).to_json(path_json, orient='records', indent=2)


json_templates = {
    'consequence_type_data.json': {
        "labels": [],
        "datasets": [{
            "label": "# of Votes",
            "data": [],
            "backgroundColor": ["#6699cc", "#fff275", "#ff8c42", "#ff3c38", "#a23e48", "#66635b", "#4e5166"],
            "borderColor": ["#6699cc", "#fff275", "#ff8c42", "#ff3c38", "#a23e48", "#66635b", "#4e5166"],
            "borderWidth": 1
        }]
    },
    'polyphen_summary_data.json': {
        "labels": [],
        "datasets": [{
            "label": "# of Votes",
            "data": [],
            "backgroundColor": ["#ff7700","#004777","#a30000","#efd28d"],
            "borderColor": ["#ff7700","#004777","#a30000","#efd28d"],
            "borderWidth": 1
        }]
    },
    'variant_class_data.json': {
        "labels": [],
        "datasets": [{
            "label": "# of Votes",
            "data": [],
            "backgroundColor": ["#ee6352","#59cd90","#3fa7d6","#fac05e","#f79d84"],
            "borderColor": ["#ee6352","#59cd90","#3fa7d6","#fac05e","#f79d84"],
            "borderWidth": 1
        }]
    }
}

labels_order = {
    'consequence_type_data.json': [ # order?
        'Intron variant',
        'Non coding transcript variant',
        'Others',
        'Intergenic variant',
        'Regulatory region variant',
        'Downstream gene variant',
        'Upstream gene variant',
    ],
}

def vep_table_to_pie_json(df, json_name):
    # Initialize JSON data structure
    json_data = json_templates[json_name].copy()
    if json_name in labels_order:
        labels = labels_order[json_name]
    else:
        df = df.sort_values('Count', ascending=False)
        labels = df['Label'].tolist()

    # Fill JSON data structure
    json_data["labels"] = []
    data_table = df.set_index('Label')['Count']
    for label in labels:
        if label not in data_table:
            continue
        count = data_table[label]
        json_data["labels"].append(f"{label} ({count:.2f}%)")
        json_data["datasets"][0]["data"].append(round(count, 2))
    return json_data


def vep_labex_tables_to_jsons(vep_dir, latex_dir, output_dir):
    if vep_dir is not None:
        vep_dir = Path(vep_dir)
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        tables = {
            'consequence_type_data.json': 'latest.consequences_table.csv',
            'polyphen_summary_data.json': 'latest.polyphen_table.csv',
            'variant_class_data.json': 'latest.var_class_table.csv',
        }
        for json_name, table_name in tables.items():
            path_table = vep_dir / table_name
            df = pd.read_csv(path_table)
            json_data = vep_table_to_pie_json(df, json_name)
            path_json = output_dir / json_name
            path_json.write_text(json.dumps(json_data, indent=2))

    if latex_dir is not None:
        latex_dir = Path(latex_dir)
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        # report_metadata.json
        output_path = output_dir / 'report_metadata.json'
        df_version = pd.read_csv(latex_dir / 'versionTable.txt', sep='\t')
        df_version['db_links'] = df_version.apply(lambda x: [{'name':x['Resource'],'link':link_warper_split(x['Website'])[1]}],axis=1)
        df_version = df_version.rename(columns={'Resource':'resource','Version':'version'})
        df_version = df_version.drop(columns=['Website'])
        df_version.to_json(output_path, orient='records', indent=2)

        # variant_calling_data.json
        output_path = output_dir / 'variant_calling_data.json'
        df_sum = pd.read_csv(latex_dir / 'latest.summary.csv',sep='\t')
        df_sum = df_sum.reset_index().rename(columns={'index':'key'})
        df_sum['key'] = df_sum['key']+1
        df_sum = df_sum.rename(columns={'Feature':'feature','Count':'count'})
        df_sum.to_json(output_path, orient='records', indent=2)

funcs = {
    'geno_tables_to_jsons': geno_tables_to_jsons,
    'vep_labex_tables_to_jsons': vep_labex_tables_to_jsons,
}

class Tools(object):
    def __init__(self):
        super(Tools, self).__init__()


if __name__ == '__main__':
    for k, v in funcs.items():
        setattr(Tools, k, staticmethod(v))
    fire.Fire(Tools)
