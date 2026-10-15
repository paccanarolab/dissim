"""OMIM 2 PubMed mapper using NCBI E-utilities instead of the OMIM API

NCBI keeps a copy of OMIM as an Entrez database, and links every OMIM entry to the
PubMed articles it cites (link name: omim_pubmed_cited). This is the same information
returned by the OMIM API referenceList endpoint, so no OMIM API key is needed.

The output has the same format as OMIM_query.py:

    OMIM ID\tPubMed ID\tPubMed ID\t...

Entries that cite no PubMed article are not written. If the output file already exists,
the MIM numbers it contains are skipped, so an interrupted run can be resumed.
"""

import configparser
import json
import time
import urllib.parse
import urllib.request

import rich.progress

ELINK_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
# requests with 200 ids are cut off by NCBI, 100 take ~10 seconds
BATCH_SIZE = 100
MAX_RETRIES = 5


def read_config(config_file):
    parser = configparser.ConfigParser()
    parser.read(config_file)
    try:
        email = parser.get("EntrezConfig", "email")
    except (configparser.NoSectionError, configparser.NoOptionError):
        raise SystemExit(f"Could not read the email from the [EntrezConfig] section of {config_file}")
    api_key = parser.get("EntrezConfig", "api_key", fallback="").strip()
    return email, api_key


def fetch_links(mim_numbers, email, api_key):
    # every MIM number must be sent as its own id parameter, otherwise
    # NCBI merges all the links in a single set and the mapping is lost
    params = {
        "dbfrom": "omim",
        "db": "pubmed",
        "linkname": "omim_pubmed_cited",
        "retmode": "json",
        "tool": "dissim",
        "email": email,
        "id": mim_numbers,
    }
    if api_key:
        params["api_key"] = api_key
    data = urllib.parse.urlencode(params, doseq=True).encode()
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            with urllib.request.urlopen(urllib.request.Request(ELINK_URL, data=data), timeout=120) as response:
                result = json.load(response)
            break
        except Exception as e:
            if attempt == MAX_RETRIES:
                raise
            print(f"Request failed ({e!r}), retrying in {5 * attempt} seconds")
            time.sleep(5 * attempt)
    if "linksets" not in result:
        raise RuntimeError(f"Unexpected response from NCBI: {result}")

    mapping = {}
    for linkset in result["linksets"]:
        pubmeds = [pmid for linksetdb in linkset.get("linksetdbs", []) for pmid in linksetdb.get("links", [])]
        if pubmeds:
            mapping[linkset["ids"][0]] = pubmeds
    return mapping


def already_done(outfile):
    try:
        with open(outfile, "r") as f:
            return {line.split("\t")[0] for line in f if line.strip()}
    except FileNotFoundError:
        return set()


def run(omim_list, outfile, config_file):
    email, api_key = read_config(config_file)
    with open(omim_list, "r") as f:
        mim_numbers = [line.strip() for line in f if line.strip()]
    done = already_done(outfile)
    if done:
        print(f"Skipping {len(done)} MIM numbers already in {outfile}")
        mim_numbers = [mim for mim in mim_numbers if mim not in done]

    # 3 requests per second without an API key, 10 with one
    delay = 0.1 if api_key else 0.34
    batches = [mim_numbers[i:i + BATCH_SIZE] for i in range(0, len(mim_numbers), BATCH_SIZE)]
    with open(outfile, "a") as o:
        for batch in rich.progress.track(batches, description="Querying NCBI..."):
            mapping = fetch_links(batch, email, api_key)
            for mim in batch:
                if mim in mapping:
                    o.write(mim + "\t" + "\t".join(mapping[mim]) + "\n")
            o.flush()
            time.sleep(delay)


if __name__ == "__main__":
    import argparse

    aparser = argparse.ArgumentParser()
    req = aparser.add_argument_group("required arguments")
    req.add_argument("--omim-list", "--in",
                     help="list of MIM numbers, one per line (output of extract_omim_list.py)",
                     required=True)
    req.add_argument("--outfile", "--o",
                     help="path to the output OMIM to PubMed mapping file",
                     required=True)
    aparser.add_argument("--config", default="PubMed2MeSH/entrez_config",
                         help="Entrez configuration file with your email and optional NCBI api_key")
    args = aparser.parse_args()
    run(args.omim_list, args.outfile, args.config)
