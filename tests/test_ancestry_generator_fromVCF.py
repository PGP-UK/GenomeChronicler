import unittest


class MyTestCase(unittest.TestCase):
    def test_something(self):
        # perl scripts/GenomeChronicler_ancestry_generator_fromVCF.pl data/NA12878wxs.vcf output/results_test_perl/ancestry_from_vcf
        # python3 scripts/GenomeChronicler_ancestry_generator_fromVCF.py data/NA12878wxs.vcf output/results_test_py/ancestry_from_vcf

        
        import os
        import sys
        from pathlib import Path
        import subprocess
        
        input_dir = Path('tests_data/ancestry_from_vcf')
        gt_dir = input_dir/'gt'
        output_dir = Path('tests_temp/ancestry_from_vcf')
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        path_gvcf = input_dir/'NA12878wxs.vcf'
        
        # os.system(f"perl perl/GenomeChronicler_ancestry_generator_fromVCF.pl {path_gvcf} {gt_dir}")
        os.system(f"python3 scripts/GenomeChronicler_ancestry_generator_fromVCF.py {path_gvcf} {output_dir}")

        files = '''
        NA12878wxs_1kGP_0-merge.missnp
        1kGP_GRCh38_nonAT_CG.bed
        1kGP_GRCh38_nonAT_CG.bed.gz
        1kGP_GRCh38_nonAT_CG.bed.gz.tbi
        #NA12878wxs.rsIDs.gvcf.gz
        NA12878wxs.nosex
        NA12878wxs.bed
        NA12878wxs.bim
        #NA12878wxs.log
        NA12878wxs.fam.mod
        NA12878wxs.fam
        NA12878wxs_1kGP_0-merge.fam
        #NA12878wxs_1kGP_0.log
        1kGP_2.nosex
        1kGP_2.bed
        1kGP_2.fam
        1kGP_2.bim
        #1kGP_2.log
        NA12878wxs_2.nosex
        NA12878wxs_2.bed
        NA12878wxs_2.fam
        NA12878wxs_2.bim
        #NA12878wxs_2.log
        NA12878wxs_1kGP.nosex
        NA12878wxs_1kGP.bed
        NA12878wxs_1kGP.fam
        NA12878wxs_1kGP.bim
        #NA12878wxs_1kGP.log
        snps_to_prune.nosex
        snps_to_prune.prune.in
        snps_to_prune.prune.out
        #snps_to_prune.log
        NA12878wxs_1kGP_pruned.nosex
        NA12878wxs_1kGP_pruned.bed
        NA12878wxs_1kGP_pruned.fam
        NA12878wxs_1kGP_pruned.bim
        #NA12878wxs_1kGP_pruned.log
        NA12878wxs_1kGP_pruned_pca_20.nosex
        NA12878wxs_1kGP_pruned_pca_20.eigenval
        NA12878wxs_1kGP_pruned_pca_20.eigenvec
        #NA12878wxs_1kGP_pruned_pca_20.log
        '''.strip().splitlines()
        for file in files:
            file = file.strip()
            gt_path = gt_dir/'temp'/file
            # if not gt_path.exists():
            if not gt_path.name.startswith('#'):
                continue
            py_path = output_dir/'temp'/file
            gt_data = gt_path.read_bytes()
            py_data = py_path.read_bytes()
            self.assertEqual(gt_data, py_data,file)

if __name__ == '__main__':
    unittest.main()
