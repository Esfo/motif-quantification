import asyncio
import time
import urllib.request
from pathlib import Path

import aiohttp
from Bio import SeqIO


folder = Path("/home/sfo/data/proteomics/fastas/")
trembend = "-NoTremb.fasta"

UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"
CRAP_URL = "ftp://ftp.thegpm.org/fasta/cRAP/crap.fasta"

uniprot_targets = [
    {
        "upo": "UP000005640",
        "name": "Human_Homo_sapien",
        "query": "proteome:UP000005640",
        "include_isoforms": True,
        "make_no_trembl": True,
    },
    {
        "upo": "UP000002254",
        "name": "Dog_Boxer_Canis_Lupus_familiaris",
        "query": "proteome:UP000002254",
        "include_isoforms": True,
        "make_no_trembl": True,
    },
    {
        "upo": "UP000002311",
        "name": "Yeast_Saccharomyces_cerevisiae",
        "query": "proteome:UP000002311",
        "include_isoforms": True,
        "make_no_trembl": True,
    },
    {
        "upo": "UP000000589",
        "name": "Mouse_Mus_musculus",
        "query": "proteome:UP000000589",
        "include_isoforms": True,
        "make_no_trembl": True,
    },
    {
        "upo": "UP000006718",
        "name": "Monkey_Rhesus_macaque_Macaca_mulatta",
        "query": "proteome:UP000006718",
        "include_isoforms": True,
        "make_no_trembl": True,
    },
    {
        # Current/reference cynomolgus macaque proteome.
        # If you specifically want the old one, change this back to UP000009130.
        "upo": "UP000233100",
        "name": "Monkey_Cynomolgus_Macaca_fascicularis",
        "query": "proteome:UP000233100",
        "include_isoforms": True,
        "make_no_trembl": True,
    },
    {
        "upo": "UP000002494",
        "name": "Rat_Rattus_norvegicus",
        "query": "proteome:UP000002494",
        "include_isoforms": True,
        "make_no_trembl": True,
    },
    {
        "upo": "UP000000625",
        "name": "Escherichia_coli",
        "query": "proteome:UP000000625",
        "include_isoforms": True,
        "make_no_trembl": True,
    },

    # Equivalent to:
    # curl -L 'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=%28reviewed%3Atrue%29%20AND%20%28organism_id%3A9606%29'
    {
        "upo": None,
        "name": "human_uniprot_reviewed",
        "query": "(reviewed:true) AND (organism_id:9606)",
        "include_isoforms": False,
        "make_no_trembl": False,
    },
]


async def fetch_uniprot_fasta(session, target):
    folder.mkdir(parents=True, exist_ok=True)

    name = target["name"]
    output_path = folder / f"{name}.fasta"
    temp_path = folder / f"{name}.fasta.tmp"

    params = {
        "compressed": "false",
        "format": "fasta",
        "query": target["query"],
    }

    if target.get("include_isoforms"):
        params["includeIsoform"] = "true"

    async with session.get(UNIPROT_STREAM_URL, params=params) as response:
        print(f"Started UniProt download for {name} - Status: {response.status}")
        response.raise_for_status()

        bytes_written = 0

        with open(temp_path, "wb") as f:
            async for chunk in response.content.iter_chunked(1024 * 1024):
                if chunk:
                    f.write(chunk)
                    bytes_written += len(chunk)

        if bytes_written == 0:
            temp_path.unlink(missing_ok=True)
            raise ValueError(f"Empty FASTA returned for {name}: {target['query']}")

        temp_path.replace(output_path)

    print(f"{name} complete - total size: {bytes_written} bytes")
    return output_path


async def download_uniprot_target(session, target, max_retries=3):
    name = target["name"]

    for attempt in range(1, max_retries + 1):
        try:
            print(f"Attempt {attempt} for {name}")
            return await fetch_uniprot_fasta(session, target)

        except asyncio.CancelledError:
            raise

        except Exception as e:
            print(f"Error downloading {name}: {e}")

            if attempt == max_retries:
                print(f"Failed to download {name} after {max_retries} attempts")
                return None

            await asyncio.sleep(5)

    return None


async def download_uniprot_files():
    timeout = aiohttp.ClientTimeout(total=600)

    async with aiohttp.ClientSession(timeout=timeout) as session:
        tasks = {
            target["name"]: asyncio.create_task(download_uniprot_target(session, target))
            for target in uniprot_targets
        }

        try:
            results = {}
            for name, task in tasks.items():
                results[name] = await task

        except asyncio.CancelledError:
            for task in tasks.values():
                task.cancel()
            await asyncio.gather(*tasks.values(), return_exceptions=True)
            raise

    failed = [name for name, path in results.items() if path is None]

    if failed:
        print(f"Failed to download from UniProt: {', '.join(failed)}")
    else:
        print("All UniProt downloads completed successfully")

    return results


def download_file_with_urllib(url, output_path):
    output_path = Path(output_path)
    temp_path = output_path.with_suffix(output_path.suffix + ".tmp")

    print(f"Started download for {output_path.name}")

    with urllib.request.urlopen(url, timeout=600) as response:
        with open(temp_path, "wb") as f:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                f.write(chunk)

    if temp_path.stat().st_size == 0:
        temp_path.unlink(missing_ok=True)
        raise ValueError(f"Empty file returned from {url}")

    temp_path.replace(output_path)

    print(f"{output_path.name} complete - total size: {output_path.stat().st_size} bytes")
    return output_path


async def download_crap_contaminants():
    output_path = folder / "contaminants_crap_gpm.fasta"
    return await asyncio.to_thread(download_file_with_urllib, CRAP_URL, output_path)


def concatenate_fastas(input_paths, output_path):
    output_path = Path(output_path)
    temp_path = output_path.with_suffix(output_path.suffix + ".tmp")

    with open(temp_path, "wb") as outfile:
        for input_path in input_paths:
            input_path = Path(input_path)

            with open(input_path, "rb") as infile:
                outfile.write(infile.read())

            outfile.write(b"\n")

    temp_path.replace(output_path)

    print(f"Wrote combined FASTA: {output_path}")
    print(f"Combined FASTA size: {output_path.stat().st_size} bytes")

    return output_path


def write_no_trembl_file(fasta_path):
    fasta_path = Path(fasta_path)
    output_path = fasta_path.with_name(f"{fasta_path.stem}{trembend}")

    reviewed_records = []

    with open(fasta_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            db = record.description.split("|", 1)[0]

            if db == "sp":
                reviewed_records.append(record)

    SeqIO.write(reviewed_records, output_path, "fasta")

    print(f"Processed {fasta_path.name}: wrote {len(reviewed_records)} reviewed Swiss-Prot records")

    return output_path


async def main():
    start_time = time.time()
    folder.mkdir(parents=True, exist_ok=True)

    print("Starting UniProt downloads...")
    uniprot_results = await download_uniprot_files()

    print("Downloading cRAP contaminants...")
    crap_path = await download_crap_contaminants()

    human_reviewed_path = uniprot_results.get("human_uniprot_reviewed")

    if human_reviewed_path is None:
        raise RuntimeError("Could not create combined FASTA because human_uniprot_reviewed.fasta failed")

    concatenate_fastas(
        [
            human_reviewed_path,
            crap_path,
        ],
        folder / "human_uniprot_reviewed_plus_contaminants.fasta",
    )

    print("Starting NoTremb file processing...")

    for target in uniprot_targets:
        if not target.get("make_no_trembl"):
            continue

        fasta_path = uniprot_results.get(target["name"])

        if fasta_path is not None:
            write_no_trembl_file(fasta_path)

    end_time = time.time()
    print(f"Total execution time: {end_time - start_time:.2f} seconds")


if __name__ == "__main__":
    asyncio.run(main())
