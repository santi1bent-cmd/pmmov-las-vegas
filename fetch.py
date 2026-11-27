"""Download the five public CDC datasets used in the analysis into data/raw/.

Sources (data.cdc.gov, Socrata API):
  j9g8-acpt  NWSS SARS-CoV-2 wastewater        site = 1179 (Flamingo WRC, Las Vegas)
  ymmh-divb  NWSS influenza A wastewater        site = 1179
  45cq-cw4i  NWSS RSV wastewater                site = 1179
  ua7e-t2fy  NHSN hospital respiratory data     jurisdiction = NV
  7mra-9cq9  NSSP emergency department visits   geography = Nevada
"""
import json, datetime, pathlib, urllib.parse, urllib.request

RAW = pathlib.Path(__file__).parent / "data" / "raw"
RAW.mkdir(parents=True, exist_ok=True)

JOBS = {
    "sars": ("j9g8-acpt", "site='1179'"),
    "flua": ("ymmh-divb", "site='1179'"),
    "rsv":  ("45cq-cw4i", "site='1179'"),
    "nhsn": ("ua7e-t2fy", "jurisdiction='NV'"),
    "nssp": ("7mra-9cq9", "geography='Nevada'"),
}

def get(dataset, where, page=50000):
    rows, offset = [], 0
    while True:
        q = urllib.parse.urlencode({"$where": where, "$limit": page, "$offset": offset, "$order": ":id"})
        url = f"https://data.cdc.gov/resource/{dataset}.json?{q}"
        with urllib.request.urlopen(url, timeout=120) as r:
            chunk = json.load(r)
        rows += chunk
        if len(chunk) < page:
            return rows
        offset += page

if __name__ == "__main__":
    for name, (ds, where) in JOBS.items():
        rows = get(ds, where)
        (RAW / f"{name}.json").write_text(json.dumps(rows))
        print(f"{name:5s} {ds}  {len(rows):6d} rows")
    (RAW / "FETCH_DATE").write_text(datetime.date.today().isoformat() + "\n")
