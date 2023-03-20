import unittest


class MyTestCase(unittest.TestCase):
    def test_something(self):
        # system("python3 ${dir}scripts/GenomeChronicler_quickFilterFinalReportTables.py ${resultsdir}/results/results_${sample}/latest.good.reportTable.csv");
        
        import os
        import sys
        from pathlib import Path

        input_dir = Path('tests_data/quick_filter_final_report_tables')
        output_dir = Path('tests_temp/quick_filter_final_report_tables')
        gt_dir = input_dir/'gt'
        Path(output_dir).mkdir(parents=True, exist_ok=True)


        import pandas as pd
        for ext in ['good','bad']:
            filename = f'latest.{ext}.reportTable.csv'
            input_path = input_dir/filename
            output_path = output_dir/filename
            output_path.parent.mkdir(parents=True,exist_ok=True)
            os.system(f"cp {input_path} {output_path}")
            os.system(f"python3 scripts/GenomeChronicler_quickFilterFinalReportTables.py {output_path}");

            gt_path = gt_dir/filename

            df_gt = pd.read_csv(gt_path)
            df_py = pd.read_csv(output_path)
            self.assertEqual(df_gt.shape,df_py.shape)
            self.assertSequenceEqual(df_gt.columns.tolist(),df_py.columns.tolist())
            for row_gt, row_py in zip(df_gt.iterrows(), df_py.iterrows()):
                row_gt = list(row_gt[1].values)
                row_py = list(row_py[1].values)
                if '\x80\x8eâ\x80\x8e' in row_gt[2]:
                    print(row_gt,row_py)
                    continue
                self.assertSequenceEqual(row_gt, row_py)

            f_gt = gt_path.read_bytes()
            f_py = output_path.read_bytes()
            self.assertEqual(f_gt, f_py, filename)

if __name__ == '__main__':
    unittest.main()
