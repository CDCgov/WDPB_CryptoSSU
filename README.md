# CryptoNet: Cryptosporidium Whole Genome Sequence Nextflow Pipeline - CryptoSSU

**Version:** 1.0.0

**Organization:** NCEZID - Division of Foodborne, Waterborne, and Environmental Diseases (DFWED) 

**Contact Email:** [ncezid_shareit@cdc.gov](mailto:ncezid_shareit@cdc.gov)  

**Description:** CryptoSSU is a **Nextflow DSL2 pipeline** designed for molecular detection and species/subtype identification of *Cryptosporidium* from whole genome (isolate) or metagenomic sequencing. Built for the CDC’s CryptoNet program, the pipeline leverages containerization (Docker/Singularity) to ensure easy deployment, reproducibility, and modular maintenance across Unix-based computing environments.

**Exemption:** Non-exempt

**Keywords:** pathogen surveillance, next generation sequencing, bioinformatics, genomics, metagenomics, parasitic disease, Cryptosporidiosis

---

> ⚠️ **Note:** The results produced by this pipeline are not ISO or CLIA-certified and should not be considered diagnostic.

---

## Overview

CryptoSSU is a **Nextflow DSL2 pipeline** designed for molecular detection and species/subtype identification of *Cryptosporidium* from whole genome (isolate) or metagenomic sequencing. Built for the CDC’s CryptoNet program, the pipeline leverages containerization (Docker/Singularity) to ensure easy deployment, reproducibility, and modular maintenance across Unix-based computing environments.

---

## Pipeline Features

The pipeline identifies *Cryptosporidium* species using a BLAST-based 18S approach and gp60 subtyping. Key processing steps include:

1. **Read Quality Control**  
   Tool: [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

2. **Adapter Trimming**  
   Tool: [`fastp`](https://github.com/OpenGene/fastp)

3. **Host/Prokaryotic Read Removal** (_optional_)  
   Tool: [`Kraken 2`](http://ccb.jhu.edu/software/kraken2/)

4. **_De novo_ Assembly**  
   Tools _(choice of one)_:  
   - [`Unicycler`](https://github.com/rrwick/Unicycler)  
   - [`SKESA`](https://github.com/ncbi/SKESA)

5. **Species Detection**
   
   Tool: Blast based approach using 18S sequence

7. **GP60 Subtype Characterization**
   
   Tool: Blast based approach using GP60 subtype sequences

8. **MultiQC Report**

   Tool: HTML or text document that summarizes the pipeline results into a report document

---

## Usage

To run the pipeline:

```bash
nextflow run CDCgov/WDPB_CryptoSSU -r main \
  -profile singularity \
  -c <config.file> \
  --outdir <output directory name> \
  --platform illumina \
  --allele_seq assets/allalleles_0.fasta \
  --input <samplesheet.csv>

```

**General disclaimer** This repository was created for use by CDC programs to collaborate on public health related projects in support of the [CDC mission](https://www.cdc.gov/about/organization/mission.htm).  GitHub is not hosted by the CDC, but is a third party website used by CDC and its partners to share information and collaborate on software. CDC use of GitHub does not imply an endorsement of any one particular service, product, or enterprise. 
  
## Public Domain Standard Notice
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](DISCLAIMER.md)
and [Code of Conduct](code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices
Please refer to [CDC's Template Repository](https://github.com/CDCgov/template) for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/main/CONTRIBUTING.md), [public domain notices and disclaimers](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md), and [code of conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md).
