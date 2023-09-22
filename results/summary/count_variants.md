# Count variants in each sample
This Python Jupyter notebook counts occurrences of each barcode in each sample from Illumina barcode sequencing, and adds these counts to the codon variant table.

## Set up analysis
### Import Python modules.
Use [plotnine](https://plotnine.readthedocs.io/en/stable/) for ggplot2-like plotting.

The analysis relies heavily on the Bloom lab's [dms_variants](https://jbloomlab.github.io/dms_variants) package:


```python
import glob
import itertools
import multiprocessing
import multiprocessing.pool
import os
import warnings

import alignparse
import alignparse.targets

import dms_variants.codonvarianttable
from dms_variants.constants import CBPALETTE
import dms_variants.illuminabarcodeparser
import dms_variants.utils
import dms_variants.plotnine_themes

from IPython.display import display, HTML

import pandas as pd

from plotnine import *

import yaml
```

Set [plotnine](https://plotnine.readthedocs.io/en/stable/) theme to the gray-grid one defined in `dms_variants`:


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using alignparse version {alignparse.__version__}")
print(f"Using dms_variants version {dms_variants.__version__}")
```

    Using alignparse version 0.6.0
    Using dms_variants version 1.4.3


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

### Parameters for notebook
Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Make output directory if needed:


```python
os.makedirs(config['counts_dir'], exist_ok=True)
os.makedirs(config['figs_dir'], exist_ok=True)
```

## Input variant tables
Initialize the table of barcode-variant pairs from the respective `process_ccs` notebooks for each background.


```python
variants = pd.read_csv(config['codon_variant_table_file_lib47'], na_filter=None)

variants = variants.reset_index(drop=True)

display(HTML(variants.tail().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>barcode</th>
      <th>variant_call_support</th>
      <th>codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>n_aa_substitutions</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>ZC45</td>
      <td>lib47_SARSr-wts</td>
      <td>TTTTTTCTAATGGAAT</td>
      <td>1</td>
      <td>NA</td>
      <td>NA</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>GD-Pangolin</td>
      <td>lib47_SARSr-wts</td>
      <td>TTTTTTCTAGCTGGAG</td>
      <td>2</td>
      <td>NA</td>
      <td>NA</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>AncSARS-CoV-1_alt</td>
      <td>lib47_SARSr-wts</td>
      <td>TTTTTTTAGAAGACAT</td>
      <td>10</td>
      <td>NA</td>
      <td>NA</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>RhGB01</td>
      <td>lib47_SARSr-wts</td>
      <td>TTTTTTTGCGTGACAT</td>
      <td>4</td>
      <td>NA</td>
      <td>NA</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>ZXC21</td>
      <td>lib47_SARSr-wts</td>
      <td>TTTTTTTGTACATAGC</td>
      <td>2</td>
      <td>NA</td>
      <td>NA</td>
      <td>0</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


Are there any barcodes in the same library that are shared across targets?
If so, we need to get rid of those as they will be confounded in barcode parsing:


```python
dup_barcodes = (
    variants
    .groupby(['library', 'barcode'])
    .size()
    .rename('duplicate_count')
    .reset_index()
    .query('duplicate_count > 1')
    )

print('Here are duplicated barcodes:')
display(HTML(dup_barcodes.head().to_html(index=False)))

print(f"\nRemoving the {len(dup_barcodes)} duplicated barcodes."
      f"Started with {len(variants)} barcodes:")
variants = (
    variants
    .merge(dup_barcodes, on=['library', 'barcode'], how='outer')
    .query('duplicate_count.isnull()', engine='python')
    )
print(f"After removing duplicates, there are {len(variants)} barcodes.")
```

    Here are duplicated barcodes:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>barcode</th>
      <th>duplicate_count</th>
    </tr>
  </thead>
  <tbody>
  </tbody>
</table>


    
    Removing the 0 duplicated barcodes.Started with 31372 barcodes:
    After removing duplicates, there are 31372 barcodes.


Pull out a target sequence for matching to the barcode and flanking sequence regions. Note, in this pipeline this is ok because our different backgrounds don't have differing flanks or other features within the actual N16 region covered in Illumina sequencing. If ever placing in-line barcodes here in the future, we would need to modify this.


```python
# get wildtype gene sequence for primary target
targets = alignparse.targets.Targets(seqsfile=config['amplicons'],
                                     feature_parse_specs=config['feature_parse_specs'])
```

## Setup to parse barcodes
Read data frame with list of all barcode runs.


```python
# barcode runs with R1 files by semicolon string split
barcode_runs = (pd.read_csv(config['barcode_runs'])
                .assign(R1=lambda x: x['R1'].str.split('; '))
                )
    
display(HTML(barcode_runs.to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>sample</th>
      <th>sample_type</th>
      <th>sort_bin</th>
      <th>concentration</th>
      <th>date</th>
      <th>number_cells</th>
      <th>R1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_01_bin1</td>
      <td>268C3</td>
      <td>1</td>
      <td>1</td>
      <td>221202</td>
      <td>474084</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s1_b1_S1_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_01_bin2</td>
      <td>268C3</td>
      <td>2</td>
      <td>1</td>
      <td>221202</td>
      <td>413453</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s1_b2_S37_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_01_bin3</td>
      <td>268C3</td>
      <td>3</td>
      <td>1</td>
      <td>221202</td>
      <td>16612</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s1_b3_S3_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_01_bin4</td>
      <td>268C3</td>
      <td>4</td>
      <td>1</td>
      <td>221202</td>
      <td>158308</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s1_b4_S38_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_02_bin1</td>
      <td>268C3</td>
      <td>1</td>
      <td>2</td>
      <td>221202</td>
      <td>653975</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s2_b1_S39_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_02_bin2</td>
      <td>268C3</td>
      <td>2</td>
      <td>2</td>
      <td>221202</td>
      <td>259584</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s2_b2_S40_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_02_bin3</td>
      <td>268C3</td>
      <td>3</td>
      <td>2</td>
      <td>221202</td>
      <td>24050</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s2_b3_S7_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_02_bin4</td>
      <td>268C3</td>
      <td>4</td>
      <td>2</td>
      <td>221202</td>
      <td>149651</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s2_b4_S8_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_03_bin1</td>
      <td>268C3</td>
      <td>1</td>
      <td>3</td>
      <td>221202</td>
      <td>816869</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s3_b1_S41_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_03_bin2</td>
      <td>268C3</td>
      <td>2</td>
      <td>3</td>
      <td>221202</td>
      <td>114509</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s3_b2_S10_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_03_bin3</td>
      <td>268C3</td>
      <td>3</td>
      <td>3</td>
      <td>221202</td>
      <td>115720</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s3_b3_S11_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_03_bin4</td>
      <td>268C3</td>
      <td>4</td>
      <td>3</td>
      <td>221202</td>
      <td>42663</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s3_b4_S12_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_04_bin1</td>
      <td>268C3</td>
      <td>1</td>
      <td>4</td>
      <td>221202</td>
      <td>899866</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s4_b1_S13_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_04_bin2</td>
      <td>268C3</td>
      <td>2</td>
      <td>4</td>
      <td>221202</td>
      <td>185960</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s4_b2_S14_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_04_bin3</td>
      <td>268C3</td>
      <td>3</td>
      <td>4</td>
      <td>221202</td>
      <td>54488</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s4_b3_S15_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_04_bin4</td>
      <td>268C3</td>
      <td>4</td>
      <td>4</td>
      <td>221202</td>
      <td>34</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s4_b4_S16_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_05_bin1</td>
      <td>268C3</td>
      <td>1</td>
      <td>5</td>
      <td>221202</td>
      <td>919590</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s5_b1_S17_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_05_bin2</td>
      <td>268C3</td>
      <td>2</td>
      <td>5</td>
      <td>221202</td>
      <td>158175</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s5_b2_S18_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_05_bin3</td>
      <td>268C3</td>
      <td>3</td>
      <td>5</td>
      <td>221202</td>
      <td>18</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s5_b3_S19_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_05_bin4</td>
      <td>268C3</td>
      <td>4</td>
      <td>5</td>
      <td>221202</td>
      <td>2</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s5_b4_S20_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_06_bin1</td>
      <td>268C3</td>
      <td>1</td>
      <td>6</td>
      <td>221202</td>
      <td>985555</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s16_b1_S61_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_06_bin2</td>
      <td>268C3</td>
      <td>2</td>
      <td>6</td>
      <td>221202</td>
      <td>82497</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s16_b2_S62_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_06_bin3</td>
      <td>268C3</td>
      <td>3</td>
      <td>6</td>
      <td>221202</td>
      <td>34</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s16_b3_S63_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C3_06_bin4</td>
      <td>268C3</td>
      <td>4</td>
      <td>6</td>
      <td>221202</td>
      <td>20</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s16_b4_S64_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_01_bin1</td>
      <td>268C183</td>
      <td>1</td>
      <td>1</td>
      <td>221202</td>
      <td>68981</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s6_b1_S21_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_01_bin2</td>
      <td>268C183</td>
      <td>2</td>
      <td>1</td>
      <td>221202</td>
      <td>99721</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s6_b2_S22_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_01_bin3</td>
      <td>268C183</td>
      <td>3</td>
      <td>1</td>
      <td>221202</td>
      <td>370462</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s6_b3_S42_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_01_bin4</td>
      <td>268C183</td>
      <td>4</td>
      <td>1</td>
      <td>221202</td>
      <td>659559</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s6_b4_S43_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_02_bin1</td>
      <td>268C183</td>
      <td>1</td>
      <td>2</td>
      <td>221202</td>
      <td>159779</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s7_b1_S25_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_02_bin2</td>
      <td>268C183</td>
      <td>2</td>
      <td>2</td>
      <td>221202</td>
      <td>216002</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s7_b2_S26_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_02_bin3</td>
      <td>268C183</td>
      <td>3</td>
      <td>2</td>
      <td>221202</td>
      <td>304450</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s7_b3_S44_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_02_bin4</td>
      <td>268C183</td>
      <td>4</td>
      <td>2</td>
      <td>221202</td>
      <td>414474</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s7_b4_S45_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_03_bin1</td>
      <td>268C183</td>
      <td>1</td>
      <td>3</td>
      <td>221202</td>
      <td>345623</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s8_b1_S46_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_03_bin2</td>
      <td>268C183</td>
      <td>2</td>
      <td>3</td>
      <td>221202</td>
      <td>294910</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s8_b2_S30_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_03_bin3</td>
      <td>268C183</td>
      <td>3</td>
      <td>3</td>
      <td>221202</td>
      <td>430819</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s8_b3_S47_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_03_bin4</td>
      <td>268C183</td>
      <td>4</td>
      <td>3</td>
      <td>221202</td>
      <td>3312</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s8_b4_S32_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_04_bin1</td>
      <td>268C183</td>
      <td>1</td>
      <td>4</td>
      <td>221202</td>
      <td>609129</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s9_b1_S48_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_04_bin2</td>
      <td>268C183</td>
      <td>2</td>
      <td>4</td>
      <td>221202</td>
      <td>452488</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s9_b2_S34_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_04_bin3</td>
      <td>268C183</td>
      <td>3</td>
      <td>4</td>
      <td>221202</td>
      <td>3129</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s9_b3_S35_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_04_bin4</td>
      <td>268C183</td>
      <td>4</td>
      <td>4</td>
      <td>221202</td>
      <td>16</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s9_b4_S36_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_05_bin1</td>
      <td>268C183</td>
      <td>1</td>
      <td>5</td>
      <td>221202</td>
      <td>905182</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s10_b1_S49_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_05_bin2</td>
      <td>268C183</td>
      <td>2</td>
      <td>5</td>
      <td>221202</td>
      <td>136770</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s10_b2_S38_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_05_bin3</td>
      <td>268C183</td>
      <td>3</td>
      <td>5</td>
      <td>221202</td>
      <td>30</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s10_b3_S39_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_05_bin4</td>
      <td>268C183</td>
      <td>4</td>
      <td>5</td>
      <td>221202</td>
      <td>13</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s10_b4_S40_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_06_bin1</td>
      <td>268C183</td>
      <td>1</td>
      <td>6</td>
      <td>221202</td>
      <td>985555</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s16_b1_S61_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_06_bin2</td>
      <td>268C183</td>
      <td>2</td>
      <td>6</td>
      <td>221202</td>
      <td>82497</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s16_b2_S62_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_06_bin3</td>
      <td>268C183</td>
      <td>3</td>
      <td>6</td>
      <td>221202</td>
      <td>34</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s16_b3_S63_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C183_06_bin4</td>
      <td>268C183</td>
      <td>4</td>
      <td>6</td>
      <td>221202</td>
      <td>20</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s16_b4_S64_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_01_bin1</td>
      <td>268C185</td>
      <td>1</td>
      <td>1</td>
      <td>221202</td>
      <td>3701</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s11_b1_S41_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_01_bin2</td>
      <td>268C185</td>
      <td>2</td>
      <td>1</td>
      <td>221202</td>
      <td>7689</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s11_b2_S42_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_01_bin3</td>
      <td>268C185</td>
      <td>3</td>
      <td>1</td>
      <td>221202</td>
      <td>111313</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s11_b3_S43_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_01_bin4</td>
      <td>268C185</td>
      <td>4</td>
      <td>1</td>
      <td>221202</td>
      <td>908426</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s11_b4_S44_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_02_bin1</td>
      <td>268C185</td>
      <td>1</td>
      <td>2</td>
      <td>221202</td>
      <td>7704</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s12_b1_S45_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_02_bin2</td>
      <td>268C185</td>
      <td>2</td>
      <td>2</td>
      <td>221202</td>
      <td>15704</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s12_b2_S46_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_02_bin3</td>
      <td>268C185</td>
      <td>3</td>
      <td>2</td>
      <td>221202</td>
      <td>122322</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s12_b3_S47_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_02_bin4</td>
      <td>268C185</td>
      <td>4</td>
      <td>2</td>
      <td>221202</td>
      <td>856538</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s12_b4_S48_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_03_bin1</td>
      <td>268C185</td>
      <td>1</td>
      <td>3</td>
      <td>221202</td>
      <td>16083</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s13_b1_S49_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_03_bin2</td>
      <td>268C185</td>
      <td>2</td>
      <td>3</td>
      <td>221202</td>
      <td>28366</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s13_b2_S50_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_03_bin3</td>
      <td>268C185</td>
      <td>3</td>
      <td>3</td>
      <td>221202</td>
      <td>669625</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s13_b3_S51_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_03_bin4</td>
      <td>268C185</td>
      <td>4</td>
      <td>3</td>
      <td>221202</td>
      <td>344125</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s13_b4_S52_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_04_bin1</td>
      <td>268C185</td>
      <td>1</td>
      <td>4</td>
      <td>221202</td>
      <td>29268</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s14_b1_S53_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_04_bin2</td>
      <td>268C185</td>
      <td>2</td>
      <td>4</td>
      <td>221202</td>
      <td>609555</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s14_b2_S50_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_04_bin3</td>
      <td>268C185</td>
      <td>3</td>
      <td>4</td>
      <td>221202</td>
      <td>401993</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s14_b3_S51_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_04_bin4</td>
      <td>268C185</td>
      <td>4</td>
      <td>4</td>
      <td>221202</td>
      <td>68</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s14_b4_S56_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_05_bin1</td>
      <td>268C185</td>
      <td>1</td>
      <td>5</td>
      <td>221202</td>
      <td>452831</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s15_b1_S57_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_05_bin2</td>
      <td>268C185</td>
      <td>2</td>
      <td>5</td>
      <td>221202</td>
      <td>562126</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s15_b2_S52_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_05_bin3</td>
      <td>268C185</td>
      <td>3</td>
      <td>5</td>
      <td>221202</td>
      <td>74</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s15_b3_S59_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_05_bin4</td>
      <td>268C185</td>
      <td>4</td>
      <td>5</td>
      <td>221202</td>
      <td>50</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s15_b4_S60_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_06_bin1</td>
      <td>268C185</td>
      <td>1</td>
      <td>6</td>
      <td>221202</td>
      <td>985555</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s16_b1_S61_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_06_bin2</td>
      <td>268C185</td>
      <td>2</td>
      <td>6</td>
      <td>221202</td>
      <td>82497</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s16_b2_S62_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_06_bin3</td>
      <td>268C185</td>
      <td>3</td>
      <td>6</td>
      <td>221202</td>
      <td>34</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s16_b3_S63_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C185_06_bin4</td>
      <td>268C185</td>
      <td>4</td>
      <td>6</td>
      <td>221202</td>
      <td>20</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s16_b4_S64_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_01_bin1</td>
      <td>268C61</td>
      <td>1</td>
      <td>1</td>
      <td>221118</td>
      <td>66324</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221118_s1_b1_S65_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_01_bin2</td>
      <td>268C61</td>
      <td>2</td>
      <td>1</td>
      <td>221118</td>
      <td>22711</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221118_s1_b2_S66_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_01_bin3</td>
      <td>268C61</td>
      <td>3</td>
      <td>1</td>
      <td>221118</td>
      <td>189973</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221118_s1_b3_S67_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_01_bin4</td>
      <td>268C61</td>
      <td>4</td>
      <td>1</td>
      <td>221118</td>
      <td>996364</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221118_s1_b4_S68_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_02_bin1</td>
      <td>268C61</td>
      <td>1</td>
      <td>2</td>
      <td>221118</td>
      <td>84182</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221118_s2_b1_S69_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_02_bin2</td>
      <td>268C61</td>
      <td>2</td>
      <td>2</td>
      <td>221118</td>
      <td>20906</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221118_s2_b2_S70_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_02_bin3</td>
      <td>268C61</td>
      <td>3</td>
      <td>2</td>
      <td>221118</td>
      <td>288821</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221118_s2_b3_S71_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_02_bin4</td>
      <td>268C61</td>
      <td>4</td>
      <td>2</td>
      <td>221118</td>
      <td>1017534</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221118_s2_b4_S72_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_03_bin1</td>
      <td>268C61</td>
      <td>1</td>
      <td>3</td>
      <td>221118</td>
      <td>80209</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221118_s3_b1_S73_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_03_bin2</td>
      <td>268C61</td>
      <td>2</td>
      <td>3</td>
      <td>221118</td>
      <td>67268</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221118_s3_b2_S74_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_03_bin3</td>
      <td>268C61</td>
      <td>3</td>
      <td>3</td>
      <td>221118</td>
      <td>738097</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221118_s3_b3_S75_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_03_bin4</td>
      <td>268C61</td>
      <td>4</td>
      <td>3</td>
      <td>221118</td>
      <td>426226</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221118_s3_b4_S76_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_04_bin1</td>
      <td>268C61</td>
      <td>1</td>
      <td>4</td>
      <td>221118</td>
      <td>112503</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221118_s4_b1_S77_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_04_bin2</td>
      <td>268C61</td>
      <td>2</td>
      <td>4</td>
      <td>221118</td>
      <td>459319</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221118_s4_b2_S78_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_04_bin3</td>
      <td>268C61</td>
      <td>3</td>
      <td>4</td>
      <td>221118</td>
      <td>725596</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221118_s4_b3_S79_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_04_bin4</td>
      <td>268C61</td>
      <td>4</td>
      <td>4</td>
      <td>221118</td>
      <td>23</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221118_s4_b4_S80_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_05_bin1</td>
      <td>268C61</td>
      <td>1</td>
      <td>5</td>
      <td>221118</td>
      <td>544792</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221118_s5_b1_S81_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_05_bin2</td>
      <td>268C61</td>
      <td>2</td>
      <td>5</td>
      <td>221118</td>
      <td>711140</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221118_s5_b2_S82_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_05_bin3</td>
      <td>268C61</td>
      <td>3</td>
      <td>5</td>
      <td>221118</td>
      <td>55</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221118_s5_b3_S83_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_05_bin4</td>
      <td>268C61</td>
      <td>4</td>
      <td>5</td>
      <td>221118</td>
      <td>13</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221118_s5_b4_S84_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_06_bin1</td>
      <td>268C61</td>
      <td>1</td>
      <td>6</td>
      <td>221118</td>
      <td>985555</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s16_b1_S61_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_06_bin2</td>
      <td>268C61</td>
      <td>2</td>
      <td>6</td>
      <td>221118</td>
      <td>82497</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s16_b2_S62_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_06_bin3</td>
      <td>268C61</td>
      <td>3</td>
      <td>6</td>
      <td>221118</td>
      <td>34</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s16_b3_S63_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>268C61_06_bin4</td>
      <td>268C61</td>
      <td>4</td>
      <td>6</td>
      <td>221118</td>
      <td>20</td>
      <td>[/uufs/chpc.utah.edu/common/home/starr-group1/sequencing/TNS/2022/221208_Overbaugh-mAb-breadth-v1/221202_s16_b4_S64_R1_001.fastq.gz]</td>
    </tr>
  </tbody>
</table>


Make sure library / sample combinations are unique:


```python
assert len(barcode_runs) == len(barcode_runs.groupby(['library', 'sample']))
```

Make sure the the libraries for which we have barcode runs are all in our variant table:


```python
unknown_libs = set(barcode_runs['library']) - set(variants['library'])
if unknown_libs:
    raise ValueError(f"Libraries with barcode runs not in variant table: {unknown_libs}")
```

Now we initialize an [IlluminaBarcodeParser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html#dms_variants.illuminabarcodeparser.IlluminaBarcodeParser) for each library.

First, get the length of the barcode from the alignment target after making sure the same length for all targets:


```python
bclen = len(targets.targets[0].get_feature('barcode').seq)

assert (bclen == len(target.get_feature('barcode').seq) for target in targets.targets)

print(f"Barcodes of length {bclen}")
```

    Barcodes of length 16


The other barcode parsing params come from the config file:


```python
parser_params = config['illumina_barcode_parser_params']

display(HTML(
    pd.Series(parser_params, name='value')
    .rename_axis(index='parameter')
    .reset_index()
    .to_html(index=False)
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>parameter</th>
      <th>value</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>upstream</td>
      <td>GGCCGC</td>
    </tr>
    <tr>
      <td>downstream</td>
      <td></td>
    </tr>
    <tr>
      <td>minq</td>
      <td>20</td>
    </tr>
    <tr>
      <td>upstream_mismatch</td>
      <td>1</td>
    </tr>
    <tr>
      <td>downstream_mismatch</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


The parser needs to know the set of valid barcodes, which are stored in the variant table and are different for each library.
So we create a different parser for each library using these valid barcode sets:


```python
# create dict keyed by library, value is parser for library
parsers = {lib: dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    bclen=bclen,
                    valid_barcodes=variants.loc[variants['library']==lib]['barcode'],
                    **parser_params)
           for lib in set(variants['library'])}

print('Number of valid barcodes searched for by each parser:')
display(HTML(
    pd.DataFrame([(lib, len(p.valid_barcodes)) for lib, p in parsers.items()],
                 columns=['library', 'number of valid barcodes'])
    .to_html(index=False)
    ))
```

    Number of valid barcodes searched for by each parser:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>number of valid barcodes</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib47_SARSr-wts</td>
      <td>31372</td>
    </tr>
  </tbody>
</table>


## Parse barcodes
We now parse the barcodes.
Since this will take a while, we utilize multiple CPUs via the Python [multiprocessing](https://docs.python.org/3.6/library/multiprocessing.html) module.
First, determine how many CPUs to use.
We use the minimum of the user-specified number hardcoded below and the number actually available.
(If you are running *interactively* on the Hutch cluster, you may need to reduce the number below in order to avoid an error as there is an enforced CPU limit on the home `rhino` nodes):


```python
ncpus = min(config['max_cpus'], multiprocessing.cpu_count())
print(f"Using {ncpus} CPUs")
```

    Using 16 CPUs


Parse the barcodes in parallel via a [multiprocessing.Pool](https://docs.python.org/3.6/library/multiprocessing.html#multiprocessing.pool.Pool) using all the available CPUs to get a list of the data frames with barcode counts / fates for each sample:


```python
def process_func(parser, r1files, library, sample):
    """Convenience function to be starmapped to multiprocessing pool."""
    return parser.parse(r1files, add_cols={'library': library, 'sample': sample})

# parallel computation of list of data frames
with multiprocessing.pool.Pool(processes=ncpus) as pool:
    bclist = pool.starmap(
                process_func,
                [(parsers[run.library], run.R1, run.library, run.sample)
                  for run in barcode_runs.itertuples()],
                )
```

Now concatenate the list into data frames of barcode counts and barcode fates:


```python
counts = pd.concat([samplecounts for samplecounts, _ in bclist],
                   sort=False,
                   ignore_index=True)

print('First few lines of counts data frame:')
display(HTML(counts.head().to_html(index=False)))

fates = pd.concat([samplefates for _, samplefates in bclist],
                  sort=False,
                  ignore_index=True)

print('First few lines of fates data frame:')
display(HTML(fates.head().to_html(index=False)))
```

    First few lines of counts data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>barcode</th>
      <th>count</th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>TATACCTACAACCAGG</td>
      <td>401</td>
      <td>lib47_SARSr-wts</td>
      <td>268C3_01_bin1</td>
    </tr>
    <tr>
      <td>AATGGAATATTCACAT</td>
      <td>392</td>
      <td>lib47_SARSr-wts</td>
      <td>268C3_01_bin1</td>
    </tr>
    <tr>
      <td>CTGCCTTACCATGTCC</td>
      <td>392</td>
      <td>lib47_SARSr-wts</td>
      <td>268C3_01_bin1</td>
    </tr>
    <tr>
      <td>GCCTACACCTAGGCTA</td>
      <td>379</td>
      <td>lib47_SARSr-wts</td>
      <td>268C3_01_bin1</td>
    </tr>
    <tr>
      <td>GCCCCCCACGAAGCTT</td>
      <td>361</td>
      <td>lib47_SARSr-wts</td>
      <td>268C3_01_bin1</td>
    </tr>
  </tbody>
</table>


    First few lines of fates data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>fate</th>
      <th>count</th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>valid barcode</td>
      <td>901532</td>
      <td>lib47_SARSr-wts</td>
      <td>268C3_01_bin1</td>
    </tr>
    <tr>
      <td>low quality barcode</td>
      <td>309435</td>
      <td>lib47_SARSr-wts</td>
      <td>268C3_01_bin1</td>
    </tr>
    <tr>
      <td>invalid barcode</td>
      <td>263039</td>
      <td>lib47_SARSr-wts</td>
      <td>268C3_01_bin1</td>
    </tr>
    <tr>
      <td>unparseable barcode</td>
      <td>25644</td>
      <td>lib47_SARSr-wts</td>
      <td>268C3_01_bin1</td>
    </tr>
    <tr>
      <td>failed chastity filter</td>
      <td>0</td>
      <td>lib47_SARSr-wts</td>
      <td>268C3_01_bin1</td>
    </tr>
  </tbody>
</table>


## Examine fates of parsed barcodes
First, we'll analyze the "fates" of the parsed barcodes.
These fates represent what happened to each Illumina read we parsed:
 - Did the barcode read fail the Illumina chastity filter?
 - Was the barcode *unparseable* (i.e., the read didn't appear to be a valid barcode based on flanking regions)?
 - Was the barcode sequence too *low quality* based on the Illumina quality scores?
 - Was the barcode parseable but *invalid* (i.e., not in our list of variant-associated barcodes in the codon variant table)?
 - Was the barcode *valid*, and so will be added to variant counts.
 
First, we just write a CSV file with all the barcode fates:


```python
fatesfile = os.path.join(config['counts_dir'], 'barcode_fates.csv')
print(f"Writing barcode fates to {fatesfile}")
fates.to_csv(fatesfile, index=False)
```

    Writing barcode fates to results/counts/barcode_fates.csv


Next, we tabulate the barcode fates in wide format:


```python
display(HTML(fates
             .pivot_table(columns='fate',
                          values='count',
                          index=['library', 'sample'])
             .to_html()
             ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>fate</th>
      <th>failed chastity filter</th>
      <th>invalid barcode</th>
      <th>low quality barcode</th>
      <th>unparseable barcode</th>
      <th>valid barcode</th>
    </tr>
    <tr>
      <th>library</th>
      <th>sample</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="96" valign="top">lib47_SARSr-wts</th>
      <th>268C183_01_bin1</th>
      <td>0</td>
      <td>72705</td>
      <td>73325</td>
      <td>5860</td>
      <td>188843</td>
    </tr>
    <tr>
      <th>268C183_01_bin2</th>
      <td>0</td>
      <td>61987</td>
      <td>74827</td>
      <td>5886</td>
      <td>205951</td>
    </tr>
    <tr>
      <th>268C183_01_bin3</th>
      <td>0</td>
      <td>496225</td>
      <td>390705</td>
      <td>25257</td>
      <td>1817477</td>
    </tr>
    <tr>
      <th>268C183_01_bin4</th>
      <td>0</td>
      <td>621878</td>
      <td>497453</td>
      <td>34788</td>
      <td>2384495</td>
    </tr>
    <tr>
      <th>268C183_02_bin1</th>
      <td>0</td>
      <td>111273</td>
      <td>127378</td>
      <td>10462</td>
      <td>350876</td>
    </tr>
    <tr>
      <th>268C183_02_bin2</th>
      <td>0</td>
      <td>56388</td>
      <td>73969</td>
      <td>6019</td>
      <td>216787</td>
    </tr>
    <tr>
      <th>268C183_02_bin3</th>
      <td>0</td>
      <td>646379</td>
      <td>481597</td>
      <td>31859</td>
      <td>2202637</td>
    </tr>
    <tr>
      <th>268C183_02_bin4</th>
      <td>0</td>
      <td>635637</td>
      <td>537211</td>
      <td>36565</td>
      <td>2591847</td>
    </tr>
    <tr>
      <th>268C183_03_bin1</th>
      <td>0</td>
      <td>616855</td>
      <td>464660</td>
      <td>31268</td>
      <td>2128751</td>
    </tr>
    <tr>
      <th>268C183_03_bin2</th>
      <td>0</td>
      <td>132832</td>
      <td>160002</td>
      <td>12763</td>
      <td>451309</td>
    </tr>
    <tr>
      <th>268C183_03_bin3</th>
      <td>0</td>
      <td>592802</td>
      <td>495479</td>
      <td>33896</td>
      <td>2371531</td>
    </tr>
    <tr>
      <th>268C183_03_bin4</th>
      <td>0</td>
      <td>4123</td>
      <td>6606</td>
      <td>536</td>
      <td>19752</td>
    </tr>
    <tr>
      <th>268C183_04_bin1</th>
      <td>0</td>
      <td>223637</td>
      <td>168403</td>
      <td>11272</td>
      <td>780107</td>
    </tr>
    <tr>
      <th>268C183_04_bin2</th>
      <td>0</td>
      <td>230816</td>
      <td>310681</td>
      <td>24810</td>
      <td>899000</td>
    </tr>
    <tr>
      <th>268C183_04_bin3</th>
      <td>0</td>
      <td>12937</td>
      <td>21444</td>
      <td>1828</td>
      <td>52505</td>
    </tr>
    <tr>
      <th>268C183_04_bin4</th>
      <td>0</td>
      <td>776</td>
      <td>618</td>
      <td>44</td>
      <td>377</td>
    </tr>
    <tr>
      <th>268C183_05_bin1</th>
      <td>0</td>
      <td>589355</td>
      <td>456843</td>
      <td>31086</td>
      <td>2121838</td>
    </tr>
    <tr>
      <th>268C183_05_bin2</th>
      <td>0</td>
      <td>81587</td>
      <td>112937</td>
      <td>9045</td>
      <td>324812</td>
    </tr>
    <tr>
      <th>268C183_05_bin3</th>
      <td>0</td>
      <td>7640</td>
      <td>2481</td>
      <td>146</td>
      <td>73</td>
    </tr>
    <tr>
      <th>268C183_05_bin4</th>
      <td>0</td>
      <td>13</td>
      <td>120</td>
      <td>1</td>
      <td>46</td>
    </tr>
    <tr>
      <th>268C183_06_bin1</th>
      <td>0</td>
      <td>612769</td>
      <td>809928</td>
      <td>64440</td>
      <td>2290784</td>
    </tr>
    <tr>
      <th>268C183_06_bin2</th>
      <td>0</td>
      <td>71760</td>
      <td>96334</td>
      <td>7177</td>
      <td>268451</td>
    </tr>
    <tr>
      <th>268C183_06_bin3</th>
      <td>0</td>
      <td>28</td>
      <td>157</td>
      <td>0</td>
      <td>69</td>
    </tr>
    <tr>
      <th>268C183_06_bin4</th>
      <td>0</td>
      <td>3</td>
      <td>3</td>
      <td>0</td>
      <td>2</td>
    </tr>
    <tr>
      <th>268C185_01_bin1</th>
      <td>0</td>
      <td>13313</td>
      <td>4818</td>
      <td>415</td>
      <td>5463</td>
    </tr>
    <tr>
      <th>268C185_01_bin2</th>
      <td>0</td>
      <td>14725</td>
      <td>11380</td>
      <td>1181</td>
      <td>28181</td>
    </tr>
    <tr>
      <th>268C185_01_bin3</th>
      <td>0</td>
      <td>85932</td>
      <td>99529</td>
      <td>7732</td>
      <td>269045</td>
    </tr>
    <tr>
      <th>268C185_01_bin4</th>
      <td>0</td>
      <td>670468</td>
      <td>909471</td>
      <td>70949</td>
      <td>2557637</td>
    </tr>
    <tr>
      <th>268C185_02_bin1</th>
      <td>0</td>
      <td>28673</td>
      <td>14318</td>
      <td>1296</td>
      <td>26446</td>
    </tr>
    <tr>
      <th>268C185_02_bin2</th>
      <td>0</td>
      <td>8453</td>
      <td>8903</td>
      <td>872</td>
      <td>23979</td>
    </tr>
    <tr>
      <th>268C185_02_bin3</th>
      <td>0</td>
      <td>92031</td>
      <td>108992</td>
      <td>8782</td>
      <td>300001</td>
    </tr>
    <tr>
      <th>268C185_02_bin4</th>
      <td>0</td>
      <td>611354</td>
      <td>817861</td>
      <td>65567</td>
      <td>2334862</td>
    </tr>
    <tr>
      <th>268C185_03_bin1</th>
      <td>0</td>
      <td>68964</td>
      <td>57252</td>
      <td>5173</td>
      <td>146226</td>
    </tr>
    <tr>
      <th>268C185_03_bin2</th>
      <td>0</td>
      <td>37585</td>
      <td>38600</td>
      <td>3167</td>
      <td>104098</td>
    </tr>
    <tr>
      <th>268C185_03_bin3</th>
      <td>0</td>
      <td>241228</td>
      <td>316731</td>
      <td>25032</td>
      <td>894201</td>
    </tr>
    <tr>
      <th>268C185_03_bin4</th>
      <td>0</td>
      <td>132676</td>
      <td>176874</td>
      <td>13944</td>
      <td>513074</td>
    </tr>
    <tr>
      <th>268C185_04_bin1</th>
      <td>0</td>
      <td>41791</td>
      <td>38193</td>
      <td>3345</td>
      <td>99126</td>
    </tr>
    <tr>
      <th>268C185_04_bin2</th>
      <td>0</td>
      <td>465247</td>
      <td>363370</td>
      <td>25322</td>
      <td>1655309</td>
    </tr>
    <tr>
      <th>268C185_04_bin3</th>
      <td>0</td>
      <td>500978</td>
      <td>402459</td>
      <td>27052</td>
      <td>1908234</td>
    </tr>
    <tr>
      <th>268C185_04_bin4</th>
      <td>0</td>
      <td>228</td>
      <td>280</td>
      <td>28</td>
      <td>751</td>
    </tr>
    <tr>
      <th>268C185_05_bin1</th>
      <td>0</td>
      <td>86137</td>
      <td>105029</td>
      <td>8277</td>
      <td>300290</td>
    </tr>
    <tr>
      <th>268C185_05_bin2</th>
      <td>0</td>
      <td>466825</td>
      <td>375666</td>
      <td>24201</td>
      <td>1751478</td>
    </tr>
    <tr>
      <th>268C185_05_bin3</th>
      <td>0</td>
      <td>39</td>
      <td>107</td>
      <td>13</td>
      <td>116</td>
    </tr>
    <tr>
      <th>268C185_05_bin4</th>
      <td>0</td>
      <td>37</td>
      <td>38</td>
      <td>1</td>
      <td>103</td>
    </tr>
    <tr>
      <th>268C185_06_bin1</th>
      <td>0</td>
      <td>612769</td>
      <td>809928</td>
      <td>64440</td>
      <td>2290784</td>
    </tr>
    <tr>
      <th>268C185_06_bin2</th>
      <td>0</td>
      <td>71760</td>
      <td>96334</td>
      <td>7177</td>
      <td>268451</td>
    </tr>
    <tr>
      <th>268C185_06_bin3</th>
      <td>0</td>
      <td>28</td>
      <td>157</td>
      <td>0</td>
      <td>69</td>
    </tr>
    <tr>
      <th>268C185_06_bin4</th>
      <td>0</td>
      <td>3</td>
      <td>3</td>
      <td>0</td>
      <td>2</td>
    </tr>
    <tr>
      <th>268C3_01_bin1</th>
      <td>0</td>
      <td>263039</td>
      <td>309435</td>
      <td>25644</td>
      <td>901532</td>
    </tr>
    <tr>
      <th>268C3_01_bin2</th>
      <td>0</td>
      <td>744397</td>
      <td>575124</td>
      <td>35834</td>
      <td>2715075</td>
    </tr>
    <tr>
      <th>268C3_01_bin3</th>
      <td>0</td>
      <td>25727</td>
      <td>38938</td>
      <td>2753</td>
      <td>108130</td>
    </tr>
    <tr>
      <th>268C3_01_bin4</th>
      <td>0</td>
      <td>479992</td>
      <td>487322</td>
      <td>34619</td>
      <td>2408136</td>
    </tr>
    <tr>
      <th>268C3_02_bin1</th>
      <td>0</td>
      <td>473990</td>
      <td>354024</td>
      <td>22454</td>
      <td>1633297</td>
    </tr>
    <tr>
      <th>268C3_02_bin2</th>
      <td>0</td>
      <td>610739</td>
      <td>477282</td>
      <td>30573</td>
      <td>2240585</td>
    </tr>
    <tr>
      <th>268C3_02_bin3</th>
      <td>0</td>
      <td>24272</td>
      <td>42034</td>
      <td>3759</td>
      <td>122478</td>
    </tr>
    <tr>
      <th>268C3_02_bin4</th>
      <td>0</td>
      <td>81310</td>
      <td>127069</td>
      <td>10513</td>
      <td>389142</td>
    </tr>
    <tr>
      <th>268C3_03_bin1</th>
      <td>0</td>
      <td>460632</td>
      <td>355244</td>
      <td>22166</td>
      <td>1621178</td>
    </tr>
    <tr>
      <th>268C3_03_bin2</th>
      <td>0</td>
      <td>73164</td>
      <td>100692</td>
      <td>8495</td>
      <td>282459</td>
    </tr>
    <tr>
      <th>268C3_03_bin3</th>
      <td>0</td>
      <td>57794</td>
      <td>89272</td>
      <td>7909</td>
      <td>261606</td>
    </tr>
    <tr>
      <th>268C3_03_bin4</th>
      <td>0</td>
      <td>17354</td>
      <td>27123</td>
      <td>2162</td>
      <td>83345</td>
    </tr>
    <tr>
      <th>268C3_04_bin1</th>
      <td>0</td>
      <td>561652</td>
      <td>708668</td>
      <td>55563</td>
      <td>1993286</td>
    </tr>
    <tr>
      <th>268C3_04_bin2</th>
      <td>0</td>
      <td>129720</td>
      <td>176057</td>
      <td>13890</td>
      <td>508851</td>
    </tr>
    <tr>
      <th>268C3_04_bin3</th>
      <td>0</td>
      <td>27170</td>
      <td>43684</td>
      <td>3595</td>
      <td>131211</td>
    </tr>
    <tr>
      <th>268C3_04_bin4</th>
      <td>0</td>
      <td>77</td>
      <td>79</td>
      <td>4</td>
      <td>208</td>
    </tr>
    <tr>
      <th>268C3_05_bin1</th>
      <td>0</td>
      <td>413466</td>
      <td>522421</td>
      <td>42040</td>
      <td>1495242</td>
    </tr>
    <tr>
      <th>268C3_05_bin2</th>
      <td>0</td>
      <td>122202</td>
      <td>169441</td>
      <td>13472</td>
      <td>485620</td>
    </tr>
    <tr>
      <th>268C3_05_bin3</th>
      <td>0</td>
      <td>26</td>
      <td>72</td>
      <td>0</td>
      <td>88</td>
    </tr>
    <tr>
      <th>268C3_05_bin4</th>
      <td>0</td>
      <td>185</td>
      <td>2060</td>
      <td>40</td>
      <td>163</td>
    </tr>
    <tr>
      <th>268C3_06_bin1</th>
      <td>0</td>
      <td>612769</td>
      <td>809928</td>
      <td>64440</td>
      <td>2290784</td>
    </tr>
    <tr>
      <th>268C3_06_bin2</th>
      <td>0</td>
      <td>71760</td>
      <td>96334</td>
      <td>7177</td>
      <td>268451</td>
    </tr>
    <tr>
      <th>268C3_06_bin3</th>
      <td>0</td>
      <td>28</td>
      <td>157</td>
      <td>0</td>
      <td>69</td>
    </tr>
    <tr>
      <th>268C3_06_bin4</th>
      <td>0</td>
      <td>3</td>
      <td>3</td>
      <td>0</td>
      <td>2</td>
    </tr>
    <tr>
      <th>268C61_01_bin1</th>
      <td>0</td>
      <td>182522</td>
      <td>254132</td>
      <td>23317</td>
      <td>735338</td>
    </tr>
    <tr>
      <th>268C61_01_bin2</th>
      <td>0</td>
      <td>48702</td>
      <td>71819</td>
      <td>5688</td>
      <td>203312</td>
    </tr>
    <tr>
      <th>268C61_01_bin3</th>
      <td>0</td>
      <td>149465</td>
      <td>166403</td>
      <td>13088</td>
      <td>458487</td>
    </tr>
    <tr>
      <th>268C61_01_bin4</th>
      <td>0</td>
      <td>784767</td>
      <td>1033856</td>
      <td>80590</td>
      <td>2882654</td>
    </tr>
    <tr>
      <th>268C61_02_bin1</th>
      <td>0</td>
      <td>69694</td>
      <td>120329</td>
      <td>9029</td>
      <td>296674</td>
    </tr>
    <tr>
      <th>268C61_02_bin2</th>
      <td>0</td>
      <td>51712</td>
      <td>57619</td>
      <td>5308</td>
      <td>160724</td>
    </tr>
    <tr>
      <th>268C61_02_bin3</th>
      <td>0</td>
      <td>252110</td>
      <td>286022</td>
      <td>22658</td>
      <td>780410</td>
    </tr>
    <tr>
      <th>268C61_02_bin4</th>
      <td>0</td>
      <td>986240</td>
      <td>1257781</td>
      <td>102348</td>
      <td>3610188</td>
    </tr>
    <tr>
      <th>268C61_03_bin1</th>
      <td>0</td>
      <td>126924</td>
      <td>181661</td>
      <td>15776</td>
      <td>534247</td>
    </tr>
    <tr>
      <th>268C61_03_bin2</th>
      <td>0</td>
      <td>139280</td>
      <td>161852</td>
      <td>11417</td>
      <td>415274</td>
    </tr>
    <tr>
      <th>268C61_03_bin3</th>
      <td>0</td>
      <td>754606</td>
      <td>913784</td>
      <td>70597</td>
      <td>2511471</td>
    </tr>
    <tr>
      <th>268C61_03_bin4</th>
      <td>0</td>
      <td>303581</td>
      <td>406238</td>
      <td>33415</td>
      <td>1186566</td>
    </tr>
    <tr>
      <th>268C61_04_bin1</th>
      <td>0</td>
      <td>64335</td>
      <td>86359</td>
      <td>7299</td>
      <td>249015</td>
    </tr>
    <tr>
      <th>268C61_04_bin2</th>
      <td>0</td>
      <td>443211</td>
      <td>509965</td>
      <td>40433</td>
      <td>1413264</td>
    </tr>
    <tr>
      <th>268C61_04_bin3</th>
      <td>0</td>
      <td>728608</td>
      <td>966319</td>
      <td>77974</td>
      <td>2741227</td>
    </tr>
    <tr>
      <th>268C61_04_bin4</th>
      <td>0</td>
      <td>13</td>
      <td>155</td>
      <td>1</td>
      <td>44</td>
    </tr>
    <tr>
      <th>268C61_05_bin1</th>
      <td>0</td>
      <td>528028</td>
      <td>633806</td>
      <td>49863</td>
      <td>1795102</td>
    </tr>
    <tr>
      <th>268C61_05_bin2</th>
      <td>0</td>
      <td>550345</td>
      <td>740362</td>
      <td>56594</td>
      <td>2039434</td>
    </tr>
    <tr>
      <th>268C61_05_bin3</th>
      <td>0</td>
      <td>43</td>
      <td>116</td>
      <td>11</td>
      <td>122</td>
    </tr>
    <tr>
      <th>268C61_05_bin4</th>
      <td>0</td>
      <td>7</td>
      <td>109</td>
      <td>0</td>
      <td>41</td>
    </tr>
    <tr>
      <th>268C61_06_bin1</th>
      <td>0</td>
      <td>612769</td>
      <td>809928</td>
      <td>64440</td>
      <td>2290784</td>
    </tr>
    <tr>
      <th>268C61_06_bin2</th>
      <td>0</td>
      <td>71760</td>
      <td>96334</td>
      <td>7177</td>
      <td>268451</td>
    </tr>
    <tr>
      <th>268C61_06_bin3</th>
      <td>0</td>
      <td>28</td>
      <td>157</td>
      <td>0</td>
      <td>69</td>
    </tr>
    <tr>
      <th>268C61_06_bin4</th>
      <td>0</td>
      <td>3</td>
      <td>3</td>
      <td>0</td>
      <td>2</td>
    </tr>
  </tbody>
</table>


Now we plot the barcode-read fates for each library / sample, showing the bars for valid barcodes in orange and the others in gray.
We see that the largest fraction of barcode reads correspond to valid barcodes, and most of the others are invalid barcodes (probably because the map to variants that aren't present in our variant table since we didn't associate all variants with barcodes). The exception to this is lib2 Titeseq_03_bin3; the PCR for this sample in the original sequencing run failed, so we followed it up with a single MiSeq lane. We did not filter out the PhiX reads from this data before parsing, so these PhiX reads will deflate the fraction of valid barcode reads as expected, but does not indicate any problems.


```python
barcode_fate_plot = (
    ggplot(
        fates
        .assign(sample=lambda x: pd.Categorical(x['sample'],
                                                x['sample'].unique(),
                                                ordered=True),
                fate=lambda x: pd.Categorical(x['fate'],
                                              x['fate'].unique(),
                                              ordered=True),
                is_valid=lambda x: x['fate'] == 'valid barcode'
                ), 
        aes('fate', 'count', fill='is_valid')) +
    geom_bar(stat='identity') +
    facet_grid('sample ~ library') +
    facet_grid('sample ~ library') +
    scale_fill_manual(CBPALETTE, guide=False) +
    theme(figure_size=(1.4 * (1 + fates['library'].nunique()),
                       1.7 * (1.2 + fates['sample'].nunique())),
          axis_text_x=element_text(angle=90),
          panel_grid_major_x=element_blank()
          ) +
    scale_y_continuous(labels=dms_variants.utils.latex_sci_not,
                       name='number of reads')
    )

_ = barcode_fate_plot.draw()
```


    
![png](count_variants_files/count_variants_42_0.png)
    


## Output csv of barcode counts in variant-barcode lookup table


```python
print(f"Writing variant counts to {config['variant_counts_file']}")
counts.to_csv(config['variant_counts_file'], index=False)
```

    Writing variant counts to results/counts/variant_counts.csv


The [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) has lots of nice functions that can be used to analyze the counts it contains.
However, we do that in the next notebook so we don't have to re-run this entire (rather computationally intensive) notebook every time we want to analyze a new aspect of the counts.


```python

```
