import unittest


class MyTestCase(unittest.TestCase):
    def test_something(self):
        import os
        import sys
        from pathlib import Path
        # os.system('')
        # self.assertEqual(True, False)  # add assertion here
        input_dir = 'tests_data/xlsx_from_tables'
        output_dir = 'tests_temp/xlsx_from_tables'
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        # os.system(f"perl perl/GenomeChronicler_XLSX_fromTables.pl ${input_dir} ${output_dir}/by_perl.xlsx")
        os.system(f"python3 scripts/GenomeChronicler_XLSX_fromTables.py {input_dir} {output_dir}/by_python.xlsx")
        path_gt = Path(input_dir)/'NA12878wxs_genotypes_23-051.xlsx'
        path_py = Path(output_dir)/'by_python.xlsx'

        import pandas as pd
        sheet_names = ['Genosets','Possibly Beneficial','Possibly Harmful']
        for sheet_name in sheet_names:
            df_gt = pd.read_excel(path_gt,sheet_name=sheet_name)
            df_py = pd.read_excel(path_py,sheet_name=sheet_name)
            self.assertEqual(df_gt.shape,df_py.shape)
            self.assertSequenceEqual(df_gt.columns.tolist(),df_py.columns.tolist())
            for row_gt, row_py in zip(df_gt.iterrows(), df_py.iterrows()):
                row_gt = list(row_gt[1].values)
                row_py = list(row_py[1].values)
                if '\x80\x8eâ\x80\x8e' in row_gt[2]:
                    print(row_gt,row_py)
                    continue
                # print(row_gt., row_py.)
                self.assertSequenceEqual(row_gt, row_py)

        # f_gt = gt_path.read_text()
        # f_py = output_path.read_text()

        # lines_gt = f_gt.splitlines()
        # lines_py = f_py.splitlines()
        # for l_gt, l_py in zip(lines_gt,lines_py):
        #     if '\u200e\u200e' in l_gt: # unicode process dfferent between perl and python
        #         continue
        #     self.assertEqual(l_gt, l_py)

        # self.assertEqual(f_perf, f_py)

if __name__ == '__main__':
    unittest.main()
