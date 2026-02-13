# DrugBank Discovery API (manual checks)

For offline verification only (API key required; not wired into DeepCoy yet):

- Base docs: https://docs.drugbank.com
- Endpoints to inspect:
  - `GET /polypeptides/<UNIPROT>/bonds` – confirm target linkage for the UniProt accession.
  - `GET /drugs/<DRUGBANK_ID>` – check `smiles` or `inchi` fields under the structure payload.
- Headers: include your API key via `Authorization: Bearer <TOKEN>` (or the scheme required by your subscription).
- When reviewing hits, prefer canonical SMILES; fall back to InChI (convert to SMILES with RDKit if needed).
- Do not scrape HTML; rely on the documented JSON endpoints only.
