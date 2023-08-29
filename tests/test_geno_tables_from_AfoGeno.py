import unittest


class MyTestCase(unittest.TestCase):
    def test_something(self):
        # system("perl ${dir}perl/GenomeChronicler_genoTables_fromAfoGeno.pl $AFOgeno_file ${resultsdir}/results/results_${sample}/ 2>>$LOGFILE2");
        # system("python3 ${dir}scripts/GenomeChronicler_genoTables_fromAfoGeno.py $AFOgeno_file ${resultsdir}/results/results_${sample}/ 2>>$LOGFILE2");

        import os
        import sys
        from pathlib import Path

        input_dir = Path('tests_data/geno_tables_from_AfoGeno')
        output_dir = Path('tests_temp/geno_tables_from_AfoGeno')
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        input_path = input_dir/'NA12878wxs.afogeno38.txt'

        # os.system(f"perl perl/GenomeChronicler_genoTables_fromAfoGeno.pl {input_path} {output_dir}/perl")
        os.system(f"python3 scripts/GenomeChronicler_genoTables_fromAfoGeno.py {input_path} {output_dir} --verbose")

        for ext in ['good','bad','genoset']:
            filename = f'latest.{ext}.reportTable.csv'

            gt_path = input_dir/'gt'/filename
            output_path = output_dir/filename

            f_gt = gt_path.read_text()
            f_py = output_path.read_text()

            lines_gt = f_gt.splitlines()
            lines_py = f_py.splitlines()
            for l_gt, l_py in zip(lines_gt,lines_py):
                if '\u200e\u200e' in l_gt: # unicode process different between perl and python
                    continue
                self.assertEqual(l_gt, l_py)

if __name__ == '__main__':
    unittest.main()
