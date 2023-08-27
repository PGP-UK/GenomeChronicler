import subprocess
from pathlib import Path

import fire


def clean_bam_file_noCHR(bam_path, output_dir=None):
    """
    Reheader the BAM file and create an index
    - Remove the "chr" prefix from the chromosome names

    Input: .bam
    Output: .clean.bam, .clean.bam.bai

    Parameters
    ----------
    bam_path: str or Path
        BAM file path
    output_dir: str, Path or None
        - If None, the output file will be saved in the same directory as the input file
        - Else, the output file will be saved in the specified directory

    Returns
    -------
    None

    """
    # Reheader the BAM file and create an index
    bam_path = Path(bam_path)
    if output_dir is None:
        output_dir = bam_path.parent
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir/ f"{bam_path.stem}.clean.BAM"
    subprocess.run(
        f"samtools reheader -c 'perl -pe \"s/^(@SQ.*)(\tSN:)chr/\$1\$2/\"' {bam_path} > {output_path}",
        shell=True)

    print(f"\t +++ INFO: Indexing the BAM file")
    subprocess.run(f"samtools index -@ 6 {output_path}", shell=True)



funcs = {
    'clean_bam_file_noCHR': clean_bam_file_noCHR,
}


class UtilsTools(object):
    def __init__(self):
        super(UtilsTools, self).__init__()


if __name__ == '__main__':
    for k, v in funcs.items():
        setattr(UtilsTools, k, staticmethod(v))
    fire.Fire(UtilsTools)


def ucfirst(s):
    if len(s) > 0:
        return s[0].upper() + s[1:]
    else:
        return s
