import unittest


class MyTestCase(unittest.TestCase):
    def test_something(self):
        # system("perl ${dir}perl/GenomeChronicler_vepTables_fromVEP.pl $VEP_file ${resultsdir}/results/results_${sample}/");

        import os
        import sys
        from pathlib import Path

        input_dir = Path('tests_data/veptable_from_vep')
        gt_dir = input_dir/'gt'
        output_dir = Path('tests_temp/veptable_from_vep')
        vep_path = input_dir/'tester_1_VEP.html'

        output_dir.mkdir(parents=True, exist_ok=True)
        os.system(f"python3 scripts/GenomeChronicler_vepTables_fromVEP.py {vep_path} {output_dir}");


        filenames = [
            'latest.chr_table.csv',
            'latest.coding_table.csv',
            'latest.consequences_table.csv',
            'latest.polyphen_table.csv',
            'latest.sift_table.csv',
            'latest.summary.csv',
            'latest.var_class_table.csv',
            'latest.var_cons_table.csv',
        ]
        for filename in filenames:
            gt_path = gt_dir/filename
            output_path = output_dir/filename

            f_gt = gt_path.read_text()
            f_py = output_path.read_text()
            self.assertEqual(f_gt, f_py)

            
if __name__ == '__main__':
    unittest.main()
