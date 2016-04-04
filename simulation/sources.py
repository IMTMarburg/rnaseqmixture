import chipseq
import ensembl
import genomics
import os
import pandas as pd
import pypipegraph as ppg
import re

pipeline = chipseq.pipelines.Illumina1_8()
aligner = chipseq.aligners.STAR()
genome = ensembl.EnsemblGenomeSQL('Homo_sapiens', 74)
fastqs_to_align = []
paired_fastqs = {}

emtab_2836_names = []
for dir in ["incoming/E-MTAB-2836"]:
    for fastq in os.listdir(dir):
        if fastq.endswith('.fastq.gz'):
            name = fastq[:fastq.rfind("_")]
            if '_1.fastq.gz' in fastq:  # we tread it as unpaired by dropping the second read...
                fastqs_to_align.append(chipseq.lanes.Lane(
                    name, os.path.join(dir, fastq), pipeline))
                emtab_2836_names.append(name)


aligned = {}
gene_exon_tag_counts = {}
for lane in fastqs_to_align:
    aligned[lane.name] = lane.align(aligner, genome)
    gene_exon_tag_counts[lane.name] = genomics.genes.annotators.ExonTagCountSmartStranded(aligned[
                                                                                          lane.name])

genes = genomics.genes.Genes(genome)
protein_coding_genes = genes.filter(
    "Protein coding", lambda df: df.gcv('biotype') == 'protein_coding')


def TAMTUMSource():
    """Our own TAMTUM (t-cell) dataset"""
    def load():
        df = pd.read_excel(
            "/martha/imt/papers/mueller/20150616_AG_Mueller_TUM_TAM_Paper/analysis_post_review/results/supplemental/RNAseq.xls",)
        df = df.set_index(['stable_id', 'name'])
        df = df[[x for x in df.columns if 'Fragment count' in x]]
        df.columns = [x[len('Fragment count TA056 '):] for x in df.columns]

        sampid_to_tissue = {}
        for x in df.columns:
            if 'TAM' in x:
                sampid_to_tissue[x] = 'TAM'
            elif 't-cell' in x.lower():
                sampid_to_tissue = 'Tcell'
            else:
                sampid_to_tissue[x] = 'TU'
        return sampid_to_tissue, df

    return 'TAMTUM', load, [ppg.FileChecksumInvariant("/martha/imt/papers/mueller/20150616_AG_Mueller_TUM_TAM_Paper/analysis_post_review/results/supplemental/RNAseq.xls")]


def GTEXSource():
    """gene tissue expression atlas dataset"""
    def load():
        # return sample_info and sample_df
        df = pd.read_csv("incoming/GTEX/GTEx_Analysis_V4_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz", sep="\t",
                         compression='gzip', skiprows=2)
        df = df.set_index(['Name'])
        df = df.drop('Description', axis=1)

        sample_df = pd.read_csv("incoming/GTEX/GTEx_Data_V4_Annotations_SampleAttributesDS.txt", sep="\t")
        samples_present = sample_df[sample_df.SAMPID.isin(df.columns)]
        sampid_to_tissue = dict((k, v) for (k, v) in zip(
            samples_present.SAMPID, samples_present.SMTSD))
        return sampid_to_tissue, df
    return 'gtex', load, []


def EMTAB_2836Source():
    def load():
        sample_df = pd.read_csv(
            "incoming/E-MTAB-2836/E-MTAB-2836.sdrf.txt", sep="\t")
        sampid_to_tissue = dict((k, v) for (k, v) in zip(
            sample_df['Comment[ENA_RUN]'], sample_df['Characteristics[organism part]']))

        df = protein_coding_genes.df.to_pandas()
        df = df[['stable_id'] +
                [gene_exon_tag_counts[x].column_name for x in emtab_2836_names]]
        df = df.set_index(['stable_id'])
        renames = {}
        for col in df.columns:
            found = re.findall("ERR[0-9]+", col)
            if found:
                renames[col] = found[0]
        df = df.rename(columns=renames)

        return sampid_to_tissue, df
    return 'emtab_2836', load, [protein_coding_genes.add_annotator(gene_exon_tag_counts[x]) for x in emtab_2836_names]

genes.write()
