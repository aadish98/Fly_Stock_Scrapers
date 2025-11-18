ortholog_mapping = {}
for result in query_results:
    query_symbol = result.get('query')
    if 'orthologs' in result:
        fly_orthologs = [h['symbol'] for h in result['homologene'] if h['taxid'] == 7227]  # 7227 = Drosophila melanogaster
        ortholog_mapping[query_symbol] = ", ".join(fly_orthologs) if fly_orthologs else None
    else:
        ortholog_mapping[query_symbol] = None  # No ortholog found

# Map to the DataFrame
df['fly_ortholog'] = df['symbol'].map(ortholog_mapping)
df.to_csv("human_genome_fly_orthologs.csv")

def human_to_fly_genes(df, symbol_col='symbol'):
    """
    Converts human gene symbols to fly gene symbols using biomaRt.

    Args:
        df: DataFrame with a column containing human gene symbols.
        symbol_col: Name of the column containing human gene symbols.

    Returns:
        DataFrame with an additional column 'fly_symbol' containing fly gene symbols.
    """

    # Connect to the Ensembl Biomart server
    server = BiomartServer("http://www.ensembl.org/biomart")
    hsapiens = server.datasets['hsapiens_gene_ensembl']
    dmelanogaster = server.datasets['dmelanogaster_gene_ensembl']

    # Define attributes to retrieve
    attributes = ['external_gene_name', 'dmelanogaster_homolog_associated_gene_name']

    # Get human gene symbols from the DataFrame
    human_genes = df[symbol_col].tolist()

    # Query the biomaRt server
    response = hsapiens.getBM(
        attributes=attributes,
        filters={'external_gene_name': human_genes},
        mart=dmelanogaster
    )

    # Create a DataFrame from the biomaRt response
    homology_df = pd.DataFrame(response, columns=['human_symbol', 'fly_symbol'])

    # Merge with the original DataFrame
    df = pd.merge(df, homology_df, left_on=symbol_col, right_on='human_symbol', how='left')

    return df