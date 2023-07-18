import unittest


class MyTestCase(unittest.TestCase):
    def test_something(self):
        # system("perl ${dir}perl/GenomeChronicler_afogeno_generator_fromBAM.pl $BAM_file $resultsdir $GATKthreads 2>>$LOGFILE2");

        
        import os
        import sys
        from pathlib import Path
        import subprocess
        
        input_dir = Path('tests_data')
        gt_dir = input_dir/'afogeno_from_bam/gt'
        output_dir = Path('tests_temp/afogeno_from_bam')
        output_dir.mkdir(parents=True, exist_ok=True)

        path_bam = input_dir/'NA12878wxs.bam.clean.BAM'
        path_bam_bai = input_dir/'NA12878wxs.bam.clean.BAM.bai' # important

        if not (gt_dir/'temp').exists():
            os.system(f"perl perl/GenomeChronicler_afogeno_generator_fromBAM.pl {path_bam} {gt_dir} 8")
            os.system(f"mv {gt_dir}/results/results_NA12878wxs/* {gt_dir}")
        # os.system(f"perl perl/GenomeChronicler_afogeno_generator_fromBAM.pl {path_bam} {output_dir} 8")
        # os.system(f"mv {output_dir}/results/results_NA12878wxs/* {output_dir}")
        os.system(f"python3 scripts/GenomeChronicler_afogeno_generator_fromBAM.py {path_bam} {output_dir} 8")
        # exit()

        # compare each line of genotypingVCF.vcf.gz
        # genotypingVCF.vcf.gz: it's different for each run
        path_gt_gvcf_gz = gt_dir/'NA12878wxs.genotypingVCF.vcf.gz'
        with subprocess.Popen(f'gzip -dcf {path_gt_gvcf_gz}', shell=True, stdout=subprocess.PIPE) as proc:
            gt_gvcf = proc.stdout.read().decode().splitlines()
        path_py_gvcf_gz = output_dir/'NA12878wxs.genotypingVCF.vcf.gz'
        with subprocess.Popen(f'gzip -dcf {path_py_gvcf_gz}', shell=True, stdout=subprocess.PIPE) as proc:
            py_gvcf = proc.stdout.read().decode().splitlines()
        self.assertEqual(len(gt_gvcf),len(py_gvcf))
        for gt_l, py_l in zip(gt_gvcf,py_gvcf):
            if 'Date=' in gt_l: # skip comparison of Date
                gt_l = gt_l.split(';')[0]
                py_l = py_l.split(';')[0]
            if 'UTC 2023' in gt_l: # skip comparison of Date
                continue
            self.assertEqual(gt_l,py_l)


        files = [
            # 'NA12878wxs.genotypes.genotypingVCF.vcf.gz',
            'temp/NA12878wxs.afogeno38.txt',
        ]        
        for file in files:
            gt_path = gt_dir/file
            py_path = output_dir/file
            gt_data = gt_path.read_bytes()
            py_data = py_path.read_bytes()
            self.assertEqual(gt_data, py_data)

if __name__ == '__main__':
    unittest.main()
