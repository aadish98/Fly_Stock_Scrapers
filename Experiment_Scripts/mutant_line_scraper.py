import requests

# Define FlyBase API endpoint for search
flybase_api_url = "https://api.flybase.org/api/v1/search"

# Define the gene symbol for Clk
gene_symbol = "Clk"

# Set up query parameters
params = {
    "query": gene_symbol,     # Search for the Clk gene
    "species": "Dmel",        # Specify Drosophila melanogaster
    "category": "alleles",    # Search for alleles
    "limit": 100              # Limit the number of results (adjust as needed)
}

# Make the request to FlyBase API
response = requests.get(flybase_api_url, params=params)

# Check if the request was successful
if response.status_code == 200:
    try:
        allele_data = response.json()  # Parse the JSON response
        print(f"Total alleles found: {allele_data['total']}")
        # Extract relevant allele and stock info
        for allele in allele_data.get('results', []):
            allele_symbol = allele.get('symbol', 'N/A')
            bdsc_stock = allele.get('stock', {}).get('BDSC', 'N/A')
            print(f"Allele: {allele_symbol}, BDSC Stock ID: {bdsc_stock}")
    except ValueError:
        print("Error parsing the response as JSON")
else:
    print(f"Error: Status Code {response.status_code}, {response.text}")
