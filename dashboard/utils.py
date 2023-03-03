link_dic = {
    "sample_accession": "https://www.ncbi.nlm.nih.gov/biosample/",
    "study_accession": "https://www.ncbi.nlm.nih.gov/bioproject/?term=",
    "run_accession": "https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=",
    "experiment_accession": "https://www.ncbi.nlm.nih.gov/sra/SRX5864477/",
    "gene": "https://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=",
    "transcript": "https://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=",
    "gene_name": "https://www.genecards.org/cgi-bin/carddisp.pl?gene="
}

def make_clickable(k: any, v: str) -> str:
    if k in link_dic:
        return f"[{v}]({link_dic[k]}{v})"
    elif v is str and v.startswith("http://"):
        return f"[{v}]({v})"
    else:
        return v

