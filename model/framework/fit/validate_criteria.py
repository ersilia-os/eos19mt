"""
One-off script to validate all ChEBI IDs used in criteria.py against the live ChEBI database.
Outputs: model/framework/fit/data/chebi_id_validation.csv

Usage: python validate_criteria.py
"""
import re
import csv
import urllib.request
import os

ROOT = os.path.dirname(os.path.abspath(__file__))
CRITERIA_FILE = os.path.join(ROOT, "..", "code", "criteria.py")
OUTPUT_FILE = os.path.join(ROOT, "data", "chebi_id_validation.csv")


def get_chebi_name(chebi_id):
    url = f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={chebi_id}"
    try:
        with urllib.request.urlopen(url, timeout=10) as resp:
            html = resp.read().decode("utf-8")
        match = re.search(r"ChEBI Name</td>\s*<td[^>]*>([^<]+)<", html)
        if match:
            return match.group(1).strip()
        match = re.search(r"<title>([^<]+)</title>", html)
        if match:
            return match.group(1).strip()
        return "NOT FOUND"
    except Exception as e:
        return f"ERROR: {e}"


def parse_criteria():
    """Return dict: chebi_id -> list of function names it appears in."""
    with open(CRITERIA_FILE) as f:
        source = f.read()

    func_blocks = re.split(r"\ndef ", source)
    result = {}
    for block in func_blocks:
        func_match = re.match(r"(\w+)\(chebis\)", block)
        if not func_match:
            continue
        func_name = func_match.group(1)
        for chebi_id in re.findall(r"CHEBI:\d+", block):
            result.setdefault(chebi_id, []).append(func_name)
    return result


if __name__ == "__main__":
    id_to_funcs = parse_criteria()
    print(f"Found {len(id_to_funcs)} unique ChEBI IDs.")

    rows = []
    for chebi_id, funcs in sorted(id_to_funcs.items(), key=lambda x: int(x[0].split(":")[1])):
        print(f"  Checking {chebi_id}...", end=" ", flush=True)
        name = get_chebi_name(chebi_id)
        print(name)
        rows.append({"chebi_id": chebi_id, "name": name, "used_in_functions": "; ".join(funcs)})

    with open(OUTPUT_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["chebi_id", "name", "used_in_functions"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nReport written to {OUTPUT_FILE}")
