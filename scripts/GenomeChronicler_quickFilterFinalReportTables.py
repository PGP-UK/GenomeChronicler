#!/usr/bin/env python3

from pathlib import Path

import pandas as pd
import fire


def filter_final_report_csv(csv_input_path, csv_output_path=None, DEBUG=False):
    '''

    '''
    # read csv
    df = pd.read_csv(csv_input_path)

    # filter
    if len(df.columns) <= 4:
        if DEBUG:
            print("WARNING: No links found")
            return
    df = df[~(df['Mag.']==0)]
    df = df[~(df['Identifier'].str.startswith('{rs')|df['Identifier'].str.startswith('{i'))]

    # export
    if csv_output_path is None: # inplace
        # warning, inplace operation
        csv_output_path = csv_input_path
    Path(csv_output_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv_output_path, index=False)


if __name__ == '__main__':
    fire.Fire(filter_final_report_csv)
