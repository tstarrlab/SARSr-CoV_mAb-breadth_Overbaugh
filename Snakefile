"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import glob
import itertools
import os.path
import os
import textwrap
import urllib.request

import Bio.SeqIO

import dms_variants.codonvarianttable
import dms_variants.illuminabarcodeparser

import pandas as pd

# Configuration  --------------------------------------------------------------
configfile: 'config.yaml'

# run "quick" rules locally:
localrules: make_dag,
            make_summary,
            save_pinned_env

# Functions -------------------------------------------------------------------
def nb_markdown(nb):
    """Return path to Markdown results of notebook `nb`."""
    return os.path.join(config['summary_dir'],
                        os.path.basename(os.path.splitext(nb)[0]) + '.md')

# Information on samples and barcode runs -------------------------------------
barcode_runs = pd.read_csv(config['barcode_runs'])

# Rules -----------------------------------------------------------------------

# making this summary is the target rule (in place of `all`) since it
# is first rule listed.
rule make_summary:
    """Create Markdown summary of analysis."""
    input:
        dag=os.path.join(config['summary_dir'], 'dag.svg'),
        env='environment_pinned.yml',
        SARSr_lib61_get_mut_bind_expr=config['SARSr_lib61_mut_bind_expr'],
        variant_counts_file=config['variant_counts_file'],
        count_variants=nb_markdown('count_variants.ipynb'),
        compute_EC50='results/summary/compute_EC50.md',
        barcode_EC50=config['mAb_EC50_file'],
        collapse_bc_lib61='results/summary/collapse_barcodes_lib61_SARSr-wts.md',
        collapse_bc_lib61_file=config['final_variant_scores_lib61_file'],
    output:
        summary = os.path.join(config['summary_dir'], 'summary.md')
    run:
        def path(f):
            """Get path relative to `summary_dir`."""
            return os.path.relpath(f, config['summary_dir'])
        with open(output.summary, 'w') as f:
            f.write(textwrap.dedent(f"""
            # Summary

            Analysis run by [Snakefile]({path(workflow.snakefile)})
            using [this config file]({path(workflow.configfiles[0])}).
            See the [README in the top directory]({path('README.md')})
            for details.

            Here is the DAG of the computational workflow:
            ![{path(input.dag)}]({path(input.dag)})

            Here is the Markdown output of each Jupyter notebook in the
            workflow:



            1. Get prior sarbecovirus homolog wildtypes expression measurements from [this repository](https://github.com/jbloomlab/SARSr-CoV_homolog_survey).

            2. [Count variants by barcode]({path(input.count_variants)}).
               Creates a [variant counts file]({path(input.variant_counts_file)})
               giving counts of each barcoded variant in each condition.

            3. [Fit EC50 to mAb binding curves]({path(input.compute_EC50)}).
               Creates a [table]({path(input.barcode_EC50)})
               giving the EC50 phenotype of each barcoded variant in each condition.
            
            4. Collapse internal replicate barcodes of each variant to final variant phenotypes for the wildtype sarbecovirus homologs pool. Analysis [here]({path(input.collapse_bc_lib61)}) and final output file [here]({path(input.collapse_bc_lib61_file)}).
                        
            """
            ).strip())

rule make_dag:
    # error message, but works: https://github.com/sequana/sequana/issues/115
    input:
        workflow.snakefile
    output:
        os.path.join(config['summary_dir'], 'dag.svg')
    shell:
        "snakemake --forceall --dag | dot -Tsvg > {output}"

rule save_pinned_env:
    input:
        dag=os.path.join(config['summary_dir'], 'dag.svg'),
    log:
    	"environment_pinned.yml"
    shell:
        """
        conda env export > {log}
        """

rule collapse_bcs_lib61_SARSr_wts:
    input:
        config['mAb_EC50_file'],
        config['SARSr_lib61_mut_bind_expr']
    output:
        config['final_variant_scores_lib61_file'],
        md='results/summary/collapse_barcodes_lib61_SARSr-wts.md',
        md_files=directory('results/summary/collapse_barcodes_lib61_SARSr-wts_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='collapse_barcodes_lib61_SARSr-wts.Rmd',
        md='collapse_barcodes_lib61_SARSr-wts.md',
        md_files='collapse_barcodes_lib61_SARSr-wts_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule calculate_bc_mAb_EC50:
    input:
        config['barcode_runs'],
        config['codon_variant_table_file_lib61'],
        config['variant_counts_file']
    output:
        config['mAb_EC50_file'],
        md='results/summary/compute_EC50.md',
        md_files=directory('results/summary/compute_EC50_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='compute_EC50.Rmd',
        md='compute_EC50.md',
        md_files='compute_EC50_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """


rule count_variants:
    """Count codon variants from Illumina barcode runs."""
    input:
        config['barcode_runs'],
        config['codon_variant_table_file_lib61']
    output:
        config['variant_counts_file'],
        nb_markdown=nb_markdown('count_variants.ipynb')
    params:
        nb='count_variants.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
rule get_SARSr_wts_mut_bind_expr:
    """Download SARSr wts library ACE2-binding and expression scores from URL."""
    output:
        file=config['SARSr_lib61_mut_bind_expr']
    run:
        urllib.request.urlretrieve(config['SARSr_lib61_mut_bind_expr_url'], output.file)