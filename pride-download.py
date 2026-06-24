#!/usr/bin/env python3

import argparse
import asyncio
import json
import os
import re
import sys
from time import monotonic
from urllib.parse import urljoin, urlparse, unquote

import aiofiles
import aiohttp
from aiohttp import ClientTimeout
from bs4 import BeautifulSoup


API_ROOT = "https://www.ebi.ac.uk/pride/ws/archive/v3"
FTP_HOST = "ftp.pride.ebi.ac.uk"
HTTPS_FTP_ROOT = f"https://{FTP_HOST}"

PROJECT_PAGES = [
    "https://www.ebi.ac.uk/pride/archive/projects/{accession}",
    "https://central.proteomexchange.org/cgi/GetDataset?ID={accession}",
    "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={accession}",
]

HEADERS = {
    "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) pridedownload/2.3",
    "Accept": "text/html,application/json,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    "Accept-Language": "en-US,en;q=0.5",
    "Connection": "keep-alive",
}

KINDS = {
    "raw": {
        "description": "Thermo RAW files",
        "file_suffixes": (".raw",),
        "dir_suffixes": (),
    },
    "thermo": {
        "description": "Thermo RAW files",
        "file_suffixes": (".raw",),
        "dir_suffixes": (),
    },
    "waters": {
        "description": "Waters .raw directories or compressed .raw uploads",
        "file_suffixes": (
            ".raw.zip",
            ".raw.tar",
            ".raw.tar.gz",
            ".raw.tgz",
            ".raw.7z",
            ".raw.rar",
        ),
        "dir_suffixes": (".raw/",),
    },
    "sciex": {
        "description": "SCIEX / AB Sciex WIFF files",
        "file_suffixes": (".wiff", ".wiff.scan", ".wiff2"),
        "dir_suffixes": (),
    },
    "bruker": {
        "description": "Bruker .d directories or compressed .d uploads",
        "file_suffixes": (
            ".d.zip",
            ".d.tar",
            ".d.tar.gz",
            ".d.tgz",
            ".d.7z",
            ".d.rar",
        ),
        "dir_suffixes": (".d/",),
    },
    "agilent": {
        "description": "Agilent .d directories or compressed .d uploads",
        "file_suffixes": (
            ".d.zip",
            ".d.tar",
            ".d.tar.gz",
            ".d.tgz",
            ".d.7z",
            ".d.rar",
        ),
        "dir_suffixes": (".d/",),
    },
    "vendor": {
        "description": "Common native vendor files/directories",
        "file_suffixes": (
            ".raw",
            ".wiff",
            ".wiff.scan",
            ".wiff2",
            ".d.zip",
            ".d.tar",
            ".d.tar.gz",
            ".d.tgz",
            ".d.7z",
            ".d.rar",
            ".raw.zip",
            ".raw.tar",
            ".raw.tar.gz",
            ".raw.tgz",
            ".raw.7z",
            ".raw.rar",
        ),
        "dir_suffixes": (".d/", ".raw/"),
    },
    "mzml": {
        "description": "mzML files",
        "file_suffixes": (".mzml", ".mzml.gz"),
        "dir_suffixes": (),
    },
    "mzxml": {
        "description": "mzXML files",
        "file_suffixes": (".mzxml", ".mzxml.gz"),
        "dir_suffixes": (),
    },
    "mgf": {
        "description": "MGF peak-list files",
        "file_suffixes": (".mgf", ".mgf.gz"),
        "dir_suffixes": (),
    },
    "mzid": {
        "description": "mzIdentML identification files",
        "file_suffixes": (".mzid", ".mzid.gz"),
        "dir_suffixes": (),
    },
    "fasta": {
        "description": "FASTA databases",
        "file_suffixes": (".fasta", ".fa", ".fas", ".faa", ".fasta.gz", ".fa.gz"),
        "dir_suffixes": (),
    },
    "txt": {
        "description": "TXT / TSV / CSV metadata or result files",
        "file_suffixes": (".txt", ".tsv", ".csv"),
        "dir_suffixes": (),
    },
    "allms": {
        "description": "Vendor raw plus mzML/mzXML/MGF",
        "file_suffixes": (
            ".raw",
            ".wiff",
            ".wiff.scan",
            ".wiff2",
            ".d.zip",
            ".d.tar",
            ".d.tar.gz",
            ".d.tgz",
            ".d.7z",
            ".d.rar",
            ".raw.zip",
            ".raw.tar",
            ".raw.tar.gz",
            ".raw.tgz",
            ".raw.7z",
            ".raw.rar",
            ".mzml",
            ".mzml.gz",
            ".mzxml",
            ".mzxml.gz",
            ".mgf",
            ".mgf.gz",
        ),
        "dir_suffixes": (".d/", ".raw/"),
    },
}


def die(message):
    raise SystemExit(message)


def normalize_accession(accession):
    accession = accession.strip().upper()

    if accession.startswith("PDX"):
        accession = "PXD" + accession[3:]

    if not re.fullmatch(r"PXD\d{6}", accession):
        die(f"expected accession like PDX######, got {accession}")

    return accession


def clean_path(url):
    return unquote(urlparse(url).path)


def lower_clean(value):
    return unquote(str(value)).split("?", 1)[0].lower()


def is_requested_file(value, kind):
    return lower_clean(value).endswith(KINDS[kind]["file_suffixes"])


def is_requested_dir(value, kind):
    value_l = lower_clean(value)
    if not value_l.endswith("/"):
        value_l += "/"
    return value_l.endswith(KINDS[kind]["dir_suffixes"])


def is_requested_target(value, kind):
    return is_requested_file(value, kind) or is_requested_dir(value, kind)


def normalize_archive_url(value, accession):
    if value is None:
        return None

    if not isinstance(value, str):
        value = json.dumps(value)

    pattern = rf"(?:https?|ftp)://{re.escape(FTP_HOST)}(/pride/data/archive/\d{{4}}/\d{{2}}/{accession}/?)"
    match = re.search(pattern, value, flags=re.I)
    if match:
        return HTTPS_FTP_ROOT + match.group(1).rstrip("/") + "/"

    pattern = rf"(/pride/data/archive/\d{{4}}/\d{{2}}/{accession}/?)"
    match = re.search(pattern, value, flags=re.I)
    if match:
        return HTTPS_FTP_ROOT + match.group(1).rstrip("/") + "/"

    return None


def normalize_file_url(value, base_url):
    value = str(value).strip().strip('"').strip("'")

    if value.startswith(f"ftp://{FTP_HOST}"):
        return value.replace(f"ftp://{FTP_HOST}", HTTPS_FTP_ROOT, 1)

    if value.startswith(f"http://{FTP_HOST}"):
        return value.replace(f"http://{FTP_HOST}", HTTPS_FTP_ROOT, 1)

    if value.startswith(f"https://{FTP_HOST}"):
        return value

    if value.startswith("/pride/data/archive/"):
        return HTTPS_FTP_ROOT + value

    return urljoin(base_url, value)


async def fetch_text(session, url, timeout):
    async with session.get(url, timeout=timeout) as response:
        response.raise_for_status()
        return await response.text()


async def fetch_json_or_text(session, url, timeout):
    text = await fetch_text(session, url, timeout)
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return text


async def resolve_archive_url(session, accession, timeout):
    api_url = f"{API_ROOT}/projects/files-path/{accession}"

    try:
        payload = await fetch_json_or_text(session, api_url, timeout)
        archive_url = normalize_archive_url(payload, accession)
        if archive_url:
            return archive_url
    except Exception as exc:
        print(f"warning: file-path API failed: {exc}", file=sys.stderr)

    for template in PROJECT_PAGES:
        page_url = template.format(accession=accession)
        try:
            html = await fetch_text(session, page_url, timeout)
            archive_url = normalize_archive_url(html, accession)
            if archive_url:
                return archive_url
        except Exception as exc:
            print(f"warning: page scrape failed for {page_url}: {exc}", file=sys.stderr)

    die(f"could not resolve PRIDE archive URL for {accession}")


def extract_file_urls_from_payload(payload, kind, archive_url):
    found = []

    def walk(obj):
        if isinstance(obj, dict):
            for key, value in obj.items():
                key_l = str(key).lower()

                if isinstance(value, str):
                    v = value.strip()

                    if FTP_HOST in v and is_requested_target(v, kind):
                        found.append(v)
                    elif v.startswith("/pride/data/archive/") and is_requested_target(v, kind):
                        found.append(v)
                    elif key_l in {"filename", "file_name", "name", "filepath", "file_path"} and is_requested_target(v, kind):
                        found.append(v)
                    elif key_l.endswith("filename") and is_requested_target(v, kind):
                        found.append(v)
                    elif key_l.endswith("filepath") and is_requested_target(v, kind):
                        found.append(v)

                walk(value)

        elif isinstance(obj, list):
            for item in obj:
                walk(item)

        elif isinstance(obj, str):
            if FTP_HOST in obj and is_requested_target(obj, kind):
                found.append(obj)

    walk(payload)

    urls = []
    seen = set()

    for value in found:
        url = normalize_file_url(value, archive_url)
        key = url.rstrip("/")
        if key not in seen:
            seen.add(key)
            urls.append(url)

    return sorted(urls, key=lambda u: clean_path(u).lower())


async def list_files_via_api(session, accession, kind, archive_url, timeout):
    filters = set(KINDS[kind]["file_suffixes"])
    filters.update(s.rstrip("/") for s in KINDS[kind]["dir_suffixes"])

    found = []

    for suffix in sorted(filters):
        api_url = f"{API_ROOT}/projects/{accession}/files/all?filenameFilter={suffix}"
        try:
            payload = await fetch_json_or_text(session, api_url, timeout)
            found.extend(extract_file_urls_from_payload(payload, kind, archive_url))
        except Exception as exc:
            print(f"warning: files API failed for filter {suffix}: {exc}", file=sys.stderr)

    return sorted(set(found), key=lambda u: clean_path(u).lower())


async def list_directory_links(session, directory_url, timeout):
    html = await fetch_text(session, directory_url, timeout)
    soup = BeautifulSoup(html, "html.parser")

    links = []
    for link in soup.find_all("a"):
        href = link.get("href")
        if not href:
            continue
        if href in {"../", "./"}:
            continue
        if href.startswith("?"):
            continue

        links.append(urljoin(directory_url, href))

    return sorted(set(links), key=lambda u: clean_path(u).lower())


async def list_top_level_matches_via_directory(session, archive_url, kind, timeout):
    links = await list_directory_links(session, archive_url, timeout)
    return [url for url in links if is_requested_target(url, kind)]


async def expand_directory_targets(session, urls, kind, timeout):
    files = []
    dirs = []

    for url in urls:
        if is_requested_dir(url, kind):
            dirs.append(url if url.endswith("/") else url + "/")
        elif is_requested_file(url, kind):
            files.append(url)

    visited_dirs = set()

    async def walk_dir(directory_url):
        directory_url = directory_url if directory_url.endswith("/") else directory_url + "/"

        if directory_url in visited_dirs:
            return

        visited_dirs.add(directory_url)

        try:
            children = await list_directory_links(session, directory_url, timeout)
        except Exception as exc:
            print(f"warning: could not list directory {directory_url}: {exc}", file=sys.stderr)
            return

        for child in children:
            if child.endswith("/"):
                await walk_dir(child)
            else:
                files.append(child)

    for directory_url in dirs:
        await walk_dir(directory_url)

    return sorted(set(files), key=lambda u: clean_path(u).lower())


def local_path_for_url(url, archive_url, project_dir):
    archive_path = clean_path(archive_url).rstrip("/") + "/"
    file_path = clean_path(url)

    if file_path.startswith(archive_path):
        relative = file_path[len(archive_path):]
    else:
        relative = os.path.basename(file_path.rstrip("/"))

    relative = relative.lstrip("/")
    return os.path.join(project_dir, relative)


async def remote_size(session, url, timeout):
    try:
        async with session.head(url, timeout=timeout, allow_redirects=True) as response:
            if response.status >= 400:
                return 0
            return int(response.headers.get("Content-Length", 0) or 0)
    except Exception:
        return 0


async def download_file(session, url, archive_url, project_dir, semaphore, timeout, chunk_size, max_retries):
    async with semaphore:
        filename = local_path_for_url(url, archive_url, project_dir)
        os.makedirs(os.path.dirname(filename) or ".", exist_ok=True)

        retries = 0

        while True:
            started = monotonic()

            try:
                local_size = os.path.getsize(filename) if os.path.exists(filename) else 0
                expected_size = await remote_size(session, url, timeout)

                if expected_size and local_size >= expected_size:
                    print(f"already finished {filename}")
                    return

                headers = {}
                mode = "ab"

                if local_size > 0:
                    headers["Range"] = f"bytes={local_size}-"

                async with session.get(url, headers=headers, timeout=timeout) as response:
                    if local_size > 0 and response.status == 200:
                        mode = "wb"
                    elif response.status not in (200, 206):
                        response.raise_for_status()

                    print(f"starting {filename}")

                    async with aiofiles.open(filename, mode) as handle:
                        async for chunk in response.content.iter_chunked(chunk_size):
                            await handle.write(chunk)

                final_size = os.path.getsize(filename) if os.path.exists(filename) else 0

                if expected_size and final_size < expected_size:
                    raise aiohttp.ClientPayloadError(
                        f"short download: {final_size} of {expected_size} bytes"
                    )

                print(f"finished {filename} {monotonic() - started:.1f}s")
                return

            except (asyncio.TimeoutError, aiohttp.ClientError, OSError) as exc:
                retries += 1

                if max_retries > 0 and retries > max_retries:
                    raise RuntimeError(f"failed {url} after {max_retries} retries") from exc

                print(f"paused {url} {exc}", file=sys.stderr)
                await asyncio.sleep(min(60, 5 * retries))
                print(f"continuing {url} retry {retries}", file=sys.stderr)


async def run(args):
    accession = normalize_accession(args.accession)
    kind = args.kind.lower()

    output_root = os.path.abspath(os.path.expanduser(args.output))
    project_dir = os.path.join(output_root, accession)
    os.makedirs(project_dir, exist_ok=True)

    timeout = ClientTimeout(
        total=args.timeout,
        sock_connect=120,
        sock_read=args.timeout,
    )

    connector = aiohttp.TCPConnector(limit=max(args.maxdownloads + 4, 16))

    async with aiohttp.ClientSession(
        headers=HEADERS,
        timeout=timeout,
        connector=connector,
    ) as session:
        archive_url = await resolve_archive_url(session, accession, timeout)
        print(f"archive {archive_url}")

        urls = await list_files_via_api(session, accession, kind, archive_url, timeout)

        if not urls:
            urls = await list_top_level_matches_via_directory(session, archive_url, kind, timeout)

        if not urls:
            die(f"no {kind} targets found for {accession}")

        download_urls = await expand_directory_targets(session, urls, kind, timeout)

        if not download_urls:
            die(f"no downloadable files found for {accession} kind={kind}")

        print(f"output {project_dir}")
        print(f"kind {kind}: {KINDS[kind]['description']}")
        print(f"found {len(download_urls)} downloadable files")

        for url in download_urls:
            print(local_path_for_url(url, archive_url, project_dir))

        if args.dry_run:
            return

        semaphore = asyncio.Semaphore(args.maxdownloads)

        await asyncio.gather(*[
            download_file(
                session=session,
                url=url,
                archive_url=archive_url,
                project_dir=project_dir,
                semaphore=semaphore,
                timeout=timeout,
                chunk_size=args.chunk_size,
                max_retries=args.max_retries,
            )
            for url in download_urls
        ])


def parse_args():
    kind_text = "\n".join(
        f"  {kind:<8} {meta['description']}"
        for kind, meta in sorted(KINDS.items())
    )

    parser = argparse.ArgumentParser(
        prog="pridedownload",
        description="Download shotgun proteomics MS files from a PRIDE project.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""kinds:
{kind_text}

examples:
  pridedownload PDX###### raw
  pridedownload PDX###### vendor
  pridedownload PDX###### mzml
  pridedownload PDX###### allms --dry-run
  pridedownload PDX###### raw -j 16
  pridedownload PDX###### raw --output /home/sfo/store/data
""",
    )

    parser.add_argument(
        "accession",
        help="PRIDE accession, e.g. PDX######. PDX###### is accepted and normalized.",
    )

    parser.add_argument(
        "kind",
        nargs="?",
        default="raw",
        choices=sorted(KINDS),
        help="file kind to download; default: raw",
    )

    parser.add_argument(
        "-o",
        "--output",
        default=".",
        help="parent output directory; default: current working directory. Files go into OUTPUT/PXD######/",
    )

    parser.add_argument(
        "-j",
        "--maxdownloads",
        type=int,
        default=8,
        help="maximum concurrent downloads; default: 8",
    )

    parser.add_argument(
        "--timeout",
        type=int,
        default=100000,
        help="aiohttp total/read timeout in seconds; default: 100000",
    )

    parser.add_argument(
        "--chunk-size",
        type=int,
        default=1024 * 1024,
        help="chunk size in bytes; default: 1048576",
    )

    parser.add_argument(
        "--max-retries",
        type=int,
        default=0,
        help="retry limit per file; 0 means retry forever",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="resolve project and print matching files without downloading",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    if args.maxdownloads < 1:
        die("--maxdownloads must be at least 1")

    asyncio.run(run(args))


if __name__ == "__main__":
    main()
