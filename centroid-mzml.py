#!/usr/bin/env python3

import argparse
import ctypes
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import threading
import traceback
import uuid
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, wait, FIRST_COMPLETED
from dataclasses import dataclass
from pathlib import Path
from time import time

import numpy as np
from pyteomics import mzml
from psims.mzml.writer import MzMLWriter


DEFAULT_DOCKER_IMAGE = "proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses:latest"
SCRIPT_VERSION = "1.9.0"

VENDOR_SUFFIXES = (
    ".raw",
    ".d",
    ".wiff",
    ".wiff2",
    ".lcd",
    ".baf",
    ".tdf",
    ".tdf_bin",
)

ACTIVE_LOCK = threading.RLock()
ACTIVE_PROCS = set()


@dataclass(frozen=True)
class PreparedInput:
    source: str
    profile: str
    centroid: str


def die(message):
    raise SystemExit(message)


def fail(message):
    raise RuntimeError(message)


def strip_ms_suffix(path):
    name = Path(path).name
    lowered = name.lower()

    for suffix in (
        ".centroid.mzml.gz",
        ".centroid.mzml",
        ".mzml.gz",
        ".mzml",
        ".raw",
        ".wiff.scan",
        ".wiff2",
        ".wiff",
        ".d",
        ".lcd",
        ".baf",
        ".tdf_bin",
        ".tdf",
    ):
        if lowered.endswith(suffix):
            return name[: -len(suffix)]

    return Path(name).stem


def is_mzml(path):
    return str(path).lower().endswith((".mzml", ".mzml.gz"))


def is_centroid_mzml(path):
    lowered = str(path).lower()
    return lowered.endswith((".centroid.mzml", ".centroid.mzml.gz"))


def is_vendor(path):
    return str(path).lower().endswith(VENDOR_SUFFIXES)


def looks_complete_mzml(path):
    path = Path(path)

    try:
        size = path.stat().st_size
    except FileNotFoundError:
        return False

    if size < 1024:
        return False

    read_size = min(1024 * 1024, size)

    with open(path, "rb") as handle:
        handle.seek(-read_size, os.SEEK_END)
        tail = handle.read().lower()

    return b"</indexedmzml>" in tail or b"</mzml>" in tail


def relative_to_or_none(path, root):
    try:
        return Path(path).resolve().relative_to(Path(root).resolve())
    except ValueError:
        return None


def choose_docker_mount_root(source, output_dir):
    source = Path(source).resolve()
    output_dir = Path(output_dir).resolve()
    cwd = Path.cwd().resolve()

    if relative_to_or_none(source, cwd) is not None and relative_to_or_none(output_dir, cwd) is not None:
        return cwd

    return Path(os.path.commonpath([str(source.parent), str(output_dir)])).resolve()


def subprocess_preexec():
    os.setsid()

    try:
        libc = ctypes.CDLL("libc.so.6", use_errno=True)
        libc.prctl(1, signal.SIGTERM)
    except Exception:
        pass

    if os.getppid() == 1:
        os.kill(os.getpid(), signal.SIGTERM)


def register_proc(proc):
    with ACTIVE_LOCK:
        ACTIVE_PROCS.add(proc)


def unregister_proc(proc):
    with ACTIVE_LOCK:
        ACTIVE_PROCS.discard(proc)


def kill_process_group(proc, timeout=1):
    if proc.poll() is not None:
        return

    if os.name == "posix":
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except ProcessLookupError:
            return
    else:
        proc.terminate()

    try:
        proc.wait(timeout=timeout)
        return
    except subprocess.TimeoutExpired:
        pass

    if os.name == "posix":
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            return
    else:
        proc.kill()

    try:
        proc.wait(timeout=timeout)
    except Exception:
        pass


def kill_active_children():
    with ACTIVE_LOCK:
        procs = list(ACTIVE_PROCS)

    for proc in procs:
        try:
            kill_process_group(proc)
        except Exception:
            pass


def hard_interrupt_exit(prepare_executor=None, centroid_executor=None):
    print("interrupted; killing active Docker/Wine/msconvert processes", file=sys.stderr, flush=True)

    try:
        kill_active_children()
    except Exception:
        pass

    if prepare_executor is not None:
        try:
            prepare_executor.shutdown(wait=False, cancel_futures=True)
        except Exception:
            pass

    if centroid_executor is not None:
        try:
            if hasattr(centroid_executor, "terminate_workers"):
                centroid_executor.terminate_workers()
            else:
                centroid_executor.shutdown(wait=False, cancel_futures=True)
        except Exception:
            pass

    os._exit(130)


def run_checked_subprocess(cmd, failure_message):
    proc = subprocess.Popen(
        cmd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        preexec_fn=subprocess_preexec if os.name == "posix" else None,
    )
    register_proc(proc)

    try:
        stdout, stderr = proc.communicate()
    except BaseException:
        kill_process_group(proc)
        raise
    finally:
        unregister_proc(proc)

    if proc.returncode != 0:
        sys.stderr.write(stdout or "")
        sys.stderr.write(stderr or "")
        fail(failure_message)

    return stdout, stderr


def run_host_msconvert(source, output_dir, args):
    source = Path(source)
    output_dir = Path(output_dir)

    cmd = [
        args.msconvert,
        str(source),
        "--mzML",
        "--64",
        "--zlib",
        f"--outdir={str(output_dir)}",
    ]

    run_checked_subprocess(cmd, f"msconvert failed for {source}")


def run_docker_msconvert(source, output_dir, args):
    source = Path(source).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    mount_root = choose_docker_mount_root(source, output_dir)
    source_rel = source.relative_to(mount_root).as_posix()
    output_rel = output_dir.relative_to(mount_root).as_posix()

    cmd = [
        "docker",
        "run",
        "--rm",
        "-e",
        "WINEDEBUG=-all",
        "-v",
        f"{mount_root}:/data",
        "-w",
        "/data",
        args.docker_image,
        *shlex.split(args.docker_msconvert),
        "--mzML",
        "--64",
        "--zlib",
        f"--outdir={output_rel}",
        source_rel,
    ]

    run_checked_subprocess(cmd, f"Docker msconvert failed for {source}")


def run_msconvert(source, output_dir, args):
    source = Path(source)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    stem = strip_ms_suffix(source)
    out_mzml = output_dir / f"{stem}.mzML"

    if out_mzml.exists() and not args.overwrite:
        if looks_complete_mzml(out_mzml):
            return out_mzml
        out_mzml.unlink()

    if args.converter == "docker":
        run_docker_msconvert(source, output_dir, args)
    elif args.converter == "host":
        run_host_msconvert(source, output_dir, args)
    else:
        fail(f"unknown converter: {args.converter}")

    if not out_mzml.exists():
        fail(f"msconvert finished but did not create {out_mzml}")

    if not looks_complete_mzml(out_mzml):
        fail(f"msconvert created incomplete mzML: {out_mzml}")

    return out_mzml


def prepare_profile_mzml(source, output_dir, args):
    source = Path(source)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    stem = strip_ms_suffix(source)
    profile_mzml = output_dir / f"{stem}.mzML"

    if is_mzml(source):
        if source.resolve() != profile_mzml.resolve():
            if args.overwrite or not profile_mzml.exists() or not looks_complete_mzml(profile_mzml):
                shutil.copy2(source, profile_mzml)
        return profile_mzml

    if is_vendor(source):
        return run_msconvert(source, output_dir, args)

    fail(f"unsupported input type: {source}")


def prepare_one(source, output_dir, args):
    source = Path(source)
    stem = strip_ms_suffix(source)
    profile_mzml = prepare_profile_mzml(source, output_dir, args)
    centroid_mzml = Path(output_dir) / f"{stem}.centroid.mzML"

    return PreparedInput(
        source=str(source),
        profile=str(profile_mzml),
        centroid=str(centroid_mzml),
    )


def derive_parallelism(args):
    if args.jobs < 1:
        die("--jobs must be at least 1")

    input_count = len(args.inputs)
    args.worker_slots = max(1, min(args.jobs, input_count))
    return args


def print_parallelism(args):
    vendor_count = sum(1 for path in args.inputs if is_vendor(path))

    if vendor_count and args.converter == "docker":
        conversion = f"up to {args.worker_slots} Docker converter(s), uncapped ProteoWizard/Wine"
    elif vendor_count:
        conversion = f"up to {args.worker_slots} host converter(s)"
    else:
        conversion = "no vendor conversion"

    print(
        f"parallelism: -j {args.jobs}; "
        f"{args.worker_slots} concurrent file job(s); "
        f"conversion: {conversion}; "
        f"centroiding: up to {args.worker_slots} worker(s)",
        file=sys.stderr,
        flush=True,
    )


def boundary_finding(fmaxes, array):
    fmaxiter = fmaxes.copy().tolist()
    fmaxiter = np.append(0, fmaxiter)
    fmaxiter = np.append(fmaxiter, len(array) - 1)

    peakbounds = []

    for n, left_anchor in enumerate(fmaxiter[:-1]):
        right_anchor = fmaxiter[n + 1] + 1

        if n > 0:
            rightseries = array[left_anchor:right_anchor]
            rightacc = np.minimum.accumulate(rightseries)
            rtrimmer = rightseries <= rightacc
            rightestimate = np.trim_zeros(rtrimmer, trim="b").size
            nr = left_anchor + rightestimate
            rightseries = array[left_anchor:nr]

            if rightseries.size:
                rcutoff = np.where(rightseries == rightseries.min())[0][0]
                rightbound = left_anchor + rcutoff + 1
            else:
                rightbound = left_anchor + 1

            peakbounds[-1].append(rightbound)

        if n < len(fmaxiter[:-1]) - 1:
            leftseries = array[left_anchor:right_anchor]
            leftacc = np.flip(np.minimum.accumulate(np.flip(leftseries)))
            ltrimmer = leftseries <= leftacc
            leftestimate = np.trim_zeros(ltrimmer, trim="f").size
            nl = right_anchor - leftestimate
            leftseries = array[nl:right_anchor]

            if leftseries.size:
                lcutoff = np.where(leftseries == leftseries.min())[0][-1]
                leftbound = nl + lcutoff
            else:
                leftbound = left_anchor

            peakbounds.append([leftbound])

    peakbounds = np.asarray(peakbounds)

    if peakbounds.size == 0:
        return np.empty((0, 3), dtype=int)

    peakparameters = np.vstack(
        (peakbounds[:, 0], fmaxes, peakbounds[:, 1])
    ).transpose()

    return np.unique(peakparameters, axis=0).astype(int)


def minpoint_reduction(barray, mindist):
    extramaxes = set()
    mask = np.repeat(False, barray.size)

    while True:
        narray = barray[~mask]

        if narray.size == 0:
            return np.array([], dtype=int)

        forwardmaxcheck = np.append(narray[:-1] > narray[1:], False)
        backwardmaxcheck = np.append(forwardmaxcheck[0], narray[1:] > narray[:-1])
        forwardmaxcheck[-1] = backwardmaxcheck[-1]

        forwardmincheck = np.append(narray[:-1] < narray[1:], False)
        backwardmincheck = np.append(forwardmincheck[0], narray[1:] < narray[:-1])
        forwardmincheck[-1] = backwardmincheck[-1]

        newmask = np.logical_and(forwardmincheck, backwardmincheck)
        mins = np.where(newmask)[0]
        maxes = np.where(np.logical_and(forwardmaxcheck, backwardmaxcheck))[0]

        extremas = np.sort(np.append(mins, maxes))

        if extremas.size == 0:
            break

        extremadistances = np.abs(extremas - extremas[:, None]) < mindist
        np.fill_diagonal(extremadistances, False)

        separatedextremas = extremas[~extremadistances.any(axis=0)]

        if separatedextremas.size > 0:
            maxestomaintain = separatedextremas[np.isin(separatedextremas, maxes)]
            maxestomaintain = (
                maxestomaintain + mask.cumsum()[~mask][maxestomaintain]
            ).tolist()
            extramaxes.update(maxestomaintain)

            minstomaintain = separatedextremas[np.isin(separatedextremas, mins)]
            newmask[minstomaintain] = False

            if minstomaintain.size > 0:
                mins = np.delete(mins, np.where(mins == minstomaintain[:, None])[1])

        adjacentextremas = extremadistances.any()

        if adjacentextremas and mins.size > 0:
            maskinds = np.argwhere(~mask)[np.argwhere(newmask)].flatten()
            mask[maskinds] = True
        else:
            break

    if not maxes.size:
        maxes = np.array([narray.argmax()])

    fmaxes = maxes + mask.cumsum()[~mask][maxes]
    fmaxes = np.unique(np.append(fmaxes, list(extramaxes))).astype(int)

    return fmaxes


def axis_peaks(array, mindist=0):
    if array.size == 0:
        return []

    maxes = minpoint_reduction(array, mindist)
    peakparameters = boundary_finding(maxes, array)

    return peakparameters.tolist()


def centroid_arrays(mza, intensityarray, intensitytype="area", masstype="average"):
    mza = np.asarray(mza, dtype=np.float64).reshape(-1)
    intensityarray = np.asarray(intensityarray, dtype=np.float64).reshape(-1)

    if mza.size == 0 or intensityarray.size == 0:
        return mza, intensityarray

    if mza.size != intensityarray.size:
        raise ValueError(f"m/z and intensity arrays differ in length: {mza.size} != {intensityarray.size}")

    peakparameters = axis_peaks(intensityarray, mindist=0)

    masses = []
    intensities = []

    for left, apex, right in peakparameters:
        pm = mza[left:right]
        pi = intensityarray[left:right]

        if pm.size == 0 or pi.size == 0:
            continue

        if masstype == "average":
            isum = pi.sum()
            if isum <= 0:
                continue
            mass = float((pm * pi).sum() / isum)
        elif masstype == "max":
            mass = float(pm[pi.argmax()])
        else:
            fail(f"bad masstype: {masstype}")

        if intensitytype == "area":
            intensity = float(np.trapezoid(pi, pm)) if pm.size > 1 else float(pi[0])
        elif intensitytype == "max":
            intensity = float(pi.max())
        elif intensitytype == "sum":
            intensity = float(pi.sum())
        else:
            fail(f"bad intensitytype: {intensitytype}")

        masses.append(mass)
        intensities.append(intensity)

    return np.asarray(masses, dtype=np.float64), np.asarray(intensities, dtype=np.float64)


def scan_ms_level(scan):
    return int(scan.get("ms level", 0) or 0)


def scan_id(scan):
    return scan.get("id") or f"index={scan.get('index', 0)}"


def scan_rt_minutes(scan):
    try:
        return float(np.real(scan["scanList"]["scan"][0]["scan start time"]))
    except Exception:
        return None


def scan_polarity(scan):
    if "negative scan" in scan:
        return "negative scan"
    if "positive scan" in scan:
        return "positive scan"
    return None


def scan_window_list(scan):
    try:
        windows = scan["scanList"]["scan"][0]["scanWindowList"]["scanWindow"]
    except Exception:
        return []

    out = []

    for window in windows:
        low = window.get("scan window lower limit")
        high = window.get("scan window upper limit")

        if low is not None and high is not None:
            out.append((float(low), float(high)))

    return out


def first_precursor_information(scan, id_map):
    try:
        precursor = scan["precursorList"]["precursor"][0]
    except Exception:
        return None

    selected = {}

    try:
        selected = precursor["selectedIonList"]["selectedIon"][0]
    except Exception:
        pass

    isolation = precursor.get("isolationWindow", {}) or {}
    activation = precursor.get("activation", {}) or {}

    info = {}

    mz_value = selected.get("selected ion m/z")
    charge = selected.get("charge state")
    intensity = selected.get("peak intensity")

    if mz_value is not None:
        info["mz"] = float(mz_value)

    if charge is not None:
        try:
            info["charge"] = int(charge)
        except Exception:
            pass

    if intensity is not None:
        info["intensity"] = float(intensity)

    spectrum_ref = precursor.get("spectrumRef")
    if spectrum_ref:
        info["spectrum_reference"] = id_map.get(spectrum_ref, spectrum_ref)

    target = isolation.get("isolation window target m/z")
    lower = isolation.get("isolation window lower offset", 0.0)
    upper = isolation.get("isolation window upper offset", 0.0)

    if target is not None:
        info["isolation_window"] = [
            float(lower or 0.0),
            float(target),
            float(upper or 0.0),
        ]

    activation_params = []

    for term in (
        "collision-induced dissociation",
        "beam-type collision-induced dissociation",
        "higher-energy collision-induced dissociation",
        "electron transfer dissociation",
        "HCD",
        "CID",
        "ETD",
    ):
        if term in activation:
            activation_params.append(term)

    if "collision energy" in activation:
        activation_params.append({"collision energy": float(activation["collision energy"])})

    if activation_params:
        info["activation"] = activation_params

    return info or None


def scan_arrays(scan):
    mz_value = scan.get("m/z array")
    intensity_value = scan.get("intensity array")

    if mz_value is None:
        mza = np.array([], dtype=np.float64)
    else:
        mza = np.asarray(mz_value, dtype=np.float64).reshape(-1)

    if intensity_value is None:
        inten = np.array([], dtype=np.float64)
    else:
        inten = np.asarray(intensity_value, dtype=np.float64).reshape(-1)

    if mza.size != inten.size:
        sid = scan_id(scan)
        raise ValueError(f"{sid}: m/z and intensity arrays differ in length: {mza.size} != {inten.size}")

    return mza, inten


def analyze_mzml(path):
    counts = defaultdict(lambda: defaultdict(int))

    with mzml.MzML(str(path), dtype=np.float64) as reader:
        for scan in reader:
            level = scan_ms_level(scan)
            counts[level]["total"] += 1

            if "profile spectrum" in scan:
                counts[level]["profile"] += 1
            elif "centroid spectrum" in scan:
                counts[level]["centroid"] += 1
            else:
                counts[level]["unknown"] += 1

    return {level: dict(stats) for level, stats in counts.items()}


def spectrum_total_from_counts(counts):
    return sum(stats.get("total", 0) for stats in counts.values())


def output_has_profile(before_counts, levels):
    for level, stats in before_counts.items():
        if level not in levels and stats.get("profile", 0):
            return True
    return False


def processing_method_for_writer(writer):
    params = [
        "MS:1000035",
        "Conversion to mzML",
    ]

    try:
        return writer.ProcessingMethod(
            order=1,
            software_reference="centroidmzml",
            params=params,
        )
    except TypeError:
        return writer.ProcessingMethod(
            order=1,
            sofware_reference="centroidmzml",
            params=params,
        )


def write_mzml_header(writer, has_profile=False):
    file_contents = [
        "MS1 spectrum",
        "MSn spectrum",
        "centroid spectrum",
    ]

    if has_profile:
        file_contents.append("profile spectrum")

    writer.file_description(file_contents=file_contents)

    writer.software_list(
        [
            {
                "id": "centroidmzml",
                "version": SCRIPT_VERSION,
                "params": ["custom unreleased software tool"],
            }
        ]
    )

    instrument = writer.InstrumentConfiguration(
        id="IC1",
        params=[],
        component_list=[],
    )
    writer.instrument_configuration_list([instrument])

    processing = writer.DataProcessing(
        [processing_method_for_writer(writer)],
        id="DP1",
    )

    writer.data_processing_list([processing])


def write_spectrum(writer, scan, spectrum_id, out_mz, out_int, out_centroided, id_map):
    level = scan_ms_level(scan)

    params = [
        {"ms level": level},
        {"total ion current": float(out_int.sum())},
    ]

    if level == 1:
        params.insert(0, "MS1 spectrum")
    elif level > 1:
        params.insert(0, "MSn spectrum")

    kwargs = {
        "mz_array": out_mz,
        "intensity_array": out_int,
        "id": spectrum_id,
        "centroided": out_centroided,
        "params": params,
        "instrument_configuration_id": "IC1",
        "compression": "zlib",
        "encoding": {
            "m/z array": np.float64,
            "intensity array": np.float64,
        },
    }

    polarity = scan_polarity(scan)
    if polarity is not None:
        kwargs["polarity"] = polarity

    precursor = first_precursor_information(scan, id_map)
    if precursor is not None:
        kwargs["precursor_information"] = precursor

    rt = scan_rt_minutes(scan)
    if rt is not None:
        kwargs["scan_start_time"] = rt

    windows = scan_window_list(scan)
    if windows:
        kwargs["scan_window_list"] = windows

    writer.write_spectrum(**kwargs)


def write_centroid_mzml(profile_mzml, centroid_mzml, levels, intensitytype, masstype, overwrite, before_counts):
    profile_mzml = Path(profile_mzml)
    centroid_mzml = Path(centroid_mzml)

    if centroid_mzml.exists() and not overwrite:
        if looks_complete_mzml(centroid_mzml):
            return
        centroid_mzml.unlink()

    if centroid_mzml.exists() and overwrite:
        centroid_mzml.unlink()

    centroid_mzml.parent.mkdir(parents=True, exist_ok=True)
    centroid_mzml.touch()

    n_spectra = spectrum_total_from_counts(before_counts)
    levels = set(levels)
    has_profile = output_has_profile(before_counts, levels)

    try:
        with MzMLWriter(open(centroid_mzml, "wb"), close=True) as writer:
            writer.controlled_vocabularies()
            write_mzml_header(writer, has_profile=has_profile)

            with writer.run(
                id=strip_ms_suffix(profile_mzml),
                instrument_configuration="IC1",
            ):
                with writer.spectrum_list(
                    count=n_spectra,
                    data_processing_method="DP1",
                ):
                    id_map = {}

                    with mzml.MzML(str(profile_mzml), dtype=np.float64) as reader:
                        for scan in reader:
                            original_id = scan_id(scan)
                            spectrum_id = original_id
                            id_map[original_id] = spectrum_id

                            level = scan_ms_level(scan)
                            mza, inten = scan_arrays(scan)

                            already_centroid = "centroid spectrum" in scan
                            unknown_centroid_state = (
                                "centroid spectrum" not in scan
                                and "profile spectrum" not in scan
                            )

                            if level in levels and not already_centroid:
                                out_mz, out_int = centroid_arrays(
                                    mza,
                                    inten,
                                    intensitytype=intensitytype,
                                    masstype=masstype,
                                )
                                out_centroided = True
                            elif level in levels and already_centroid:
                                out_mz, out_int = mza, inten
                                out_centroided = True
                            else:
                                out_mz, out_int = mza, inten
                                out_centroided = already_centroid and not unknown_centroid_state

                            write_spectrum(
                                writer=writer,
                                scan=scan,
                                spectrum_id=spectrum_id,
                                out_mz=out_mz,
                                out_int=out_int,
                                out_centroided=out_centroided,
                                id_map=id_map,
                            )

    except Exception:
        try:
            centroid_mzml.unlink()
        except FileNotFoundError:
            pass
        raise


def classification_for_level(stats):
    total = stats.get("total", 0)
    profile = stats.get("profile", 0)
    centroid = stats.get("centroid", 0)
    unknown = stats.get("unknown", 0)

    if total == 0:
        return "absent"
    if profile == total:
        return "profile"
    if centroid == total:
        return "centroid"
    if unknown == total:
        return "unknown"

    return "mixed"


def make_project_notes(source, profile_mzml, centroid_mzml, before_counts, levels, converter):
    lines = []
    lines.append("=" * 80)
    lines.append(f"source: {source}")
    lines.append(f"converter: {converter}")
    lines.append(f"profile mzML: {profile_mzml}")
    lines.append(f"centroid mzML: {centroid_mzml}")
    lines.append("spectrum IDs in centroid mzML: preserved from profile mzML")
    lines.append("")

    explicit_profile = False

    for level in sorted(before_counts):
        stats = before_counts[level]
        explicit_profile = explicit_profile or stats.get("profile", 0) > 0

        lines.append(
            "MS{level}: total={total} profile={profile} centroid={centroid} unknown={unknown} classification={kind}".format(
                level=level,
                total=stats.get("total", 0),
                profile=stats.get("profile", 0),
                centroid=stats.get("centroid", 0),
                unknown=stats.get("unknown", 0),
                kind=classification_for_level(stats),
            )
        )

    lines.append("")

    for level in sorted(levels):
        stats = before_counts.get(level, {})
        kind = classification_for_level(stats)

        if kind == "profile":
            lines.append(f"MS{level}: profile input; custom centroiding applied.")
        elif kind == "centroid":
            lines.append(f"MS{level}: already centroided input; copied into .centroid.mzML.")
        elif kind == "mixed":
            lines.append(f"MS{level}: mixed profile/centroid input; non-centroid scans centroided, centroid scans copied.")
        elif kind == "unknown":
            lines.append(f"MS{level}: mzML did not explicitly say profile or centroid; treated as needing centroiding.")
        elif kind == "absent":
            lines.append(f"MS{level}: no scans found.")

    if not explicit_profile:
        lines.append("")
        lines.append("No scans were explicitly marked as profile spectrum. The profile mzML may already contain centroided data.")

    lines.append("")
    return "\n".join(lines) + "\n"


def process_prepared(prepared, levels, intensity, mass, overwrite, converter):
    before_counts = analyze_mzml(prepared.profile)

    write_centroid_mzml(
        profile_mzml=prepared.profile,
        centroid_mzml=prepared.centroid,
        levels=set(levels),
        intensitytype=intensity,
        masstype=mass,
        overwrite=overwrite,
        before_counts=before_counts,
    )

    notes = make_project_notes(
        source=prepared.source,
        profile_mzml=prepared.profile,
        centroid_mzml=prepared.centroid,
        before_counts=before_counts,
        levels=set(levels),
        converter=converter,
    )

    return {
        "source": prepared.source,
        "profile": prepared.profile,
        "centroid": prepared.centroid,
        "notes": notes,
    }


def parse_levels(value):
    levels = []

    for part in value.split(","):
        part = part.strip()

        if not part:
            continue

        if "-" in part:
            left, right = part.split("-", 1)
            levels.extend(range(int(left), int(right) + 1))
        else:
            levels.append(int(part))

    if not levels:
        die("empty --levels")

    return sorted(set(levels))


def parse_args():
    default_jobs = max(1, min(4, os.cpu_count() or 1))

    parser = argparse.ArgumentParser(
        prog="centroid-mzml",
        description="Convert vendor/profile MS data to profile mzML and custom-centroided mzML.",
    )

    parser.add_argument("inputs", nargs="+", help="vendor files/directories or mzML files")

    parser.add_argument(
        "-o",
        "--output",
        default=".",
        help="output directory; default: current directory",
    )

    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=default_jobs,
        help=f"concurrent file jobs; default: {default_jobs}",
    )

    parser.add_argument(
        "--levels",
        type=parse_levels,
        default=parse_levels("1,2"),
        help="MS levels to centroid; default: 1,2",
    )

    parser.add_argument(
        "--intensity",
        choices=("area", "max", "sum"),
        default="area",
        help="centroid intensity calculation; default: area",
    )

    parser.add_argument(
        "--mass",
        choices=("average", "max"),
        default="average",
        help="centroid m/z calculation; default: average",
    )

    parser.add_argument(
        "--converter",
        choices=("docker", "host"),
        default="docker",
        help="vendor converter backend; default: docker",
    )

    parser.add_argument(
        "--docker-image",
        default=DEFAULT_DOCKER_IMAGE,
        help=f"Docker ProteoWizard image; default: {DEFAULT_DOCKER_IMAGE}",
    )

    parser.add_argument(
        "--docker-msconvert",
        default="wine msconvert",
        help='command inside Docker image; default: "wine msconvert"',
    )

    parser.add_argument(
        "--msconvert",
        default="msconvert",
        help="host msconvert executable, only used with --converter host",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="overwrite existing .mzML and .centroid.mzML outputs",
    )

    parser.add_argument(
        "--tracebacks",
        action="store_true",
        help="print full Python tracebacks for per-file failures",
    )

    return parser.parse_args()


def filter_inputs(inputs):
    kept = []
    skipped = 0

    for path in inputs:
        if is_centroid_mzml(path):
            skipped += 1
            continue
        kept.append(path)

    if skipped:
        print(f"skipped {skipped} existing centroid mzML input(s)", file=sys.stderr, flush=True)

    if not kept:
        die("no non-centroid inputs left to process")

    return kept


def write_notes(output_dir, notes_blocks):
    notes_path = Path(output_dir) / "project-notes.txt"

    with open(notes_path, "a", encoding="utf-8") as handle:
        for block in notes_blocks:
            handle.write(block)


def validate_inputs(inputs):
    missing = [str(path) for path in inputs if not Path(path).exists()]

    if missing:
        die("missing input(s): " + ", ".join(missing))


def submit_centroid(centroid_executor, prepared, args):
    centroid_path = Path(prepared.centroid)

    if not centroid_path.exists():
        centroid_path.parent.mkdir(parents=True, exist_ok=True)
        centroid_path.touch()

    print(f"centroiding {prepared.centroid}", file=sys.stderr, flush=True)

    return centroid_executor.submit(
        process_prepared,
        prepared,
        args.levels,
        args.intensity,
        args.mass,
        args.overwrite,
        args.converter,
    )


def print_failure(label, item, exc, args):
    print(f"failed {label} {item}: {exc}", file=sys.stderr, flush=True)

    if args.tracebacks:
        print(traceback.format_exc(), file=sys.stderr)


def run_pipeline(args, output_dir):
    validate_inputs(args.inputs)

    results = []
    notes_blocks = []
    failures = []

    prepare_executor = None
    centroid_executor = None

    pending_sources = list(args.inputs)
    prepare_futures = {}
    centroid_futures = {}

    try:
        prepare_executor = ThreadPoolExecutor(max_workers=args.worker_slots)
        centroid_executor = ProcessPoolExecutor(max_workers=args.worker_slots)

        while pending_sources or prepare_futures or centroid_futures:
            active_jobs = len(prepare_futures) + len(centroid_futures)
            free_slots = args.worker_slots - active_jobs

            while free_slots > 0 and pending_sources:
                source = pending_sources.pop(0)
                future = prepare_executor.submit(prepare_one, source, output_dir, args)
                prepare_futures[future] = source
                free_slots -= 1

            watched = set(prepare_futures) | set(centroid_futures)

            if not watched:
                continue

            done, _ = wait(watched, timeout=0.2, return_when=FIRST_COMPLETED)

            if not done:
                continue

            for future in done:
                if future in prepare_futures:
                    source = prepare_futures.pop(future)

                    try:
                        prepared = future.result()
                    except Exception as exc:
                        failures.append((source, exc))
                        print_failure("preparing", source, exc, args)
                        continue

                    print(f"profile   {prepared.profile}", flush=True)

                    try:
                        centroid_future = submit_centroid(centroid_executor, prepared, args)
                    except Exception as exc:
                        failures.append((prepared.source, exc))
                        print_failure("submitting centroid job for", prepared.source, exc, args)
                        continue

                    centroid_futures[centroid_future] = prepared
                    continue

                if future in centroid_futures:
                    prepared = centroid_futures.pop(future)

                    try:
                        result = future.result()
                    except Exception as exc:
                        failures.append((prepared.source, exc))
                        print_failure("centroiding", prepared.source, exc, args)
                        continue

                    results.append(result)
                    notes_blocks.append(result["notes"])
                    print(f"centroided {result['centroid']}", flush=True)

    except KeyboardInterrupt:
        hard_interrupt_exit(prepare_executor, centroid_executor)

    finally:
        if prepare_executor is not None:
            prepare_executor.shutdown(wait=False, cancel_futures=True)
        if centroid_executor is not None:
            centroid_executor.shutdown(wait=False, cancel_futures=True)

    if notes_blocks:
        write_notes(output_dir, notes_blocks)

    if failures:
        print(
            f"finished with {len(failures)} failed file(s) and {len(results)} successful centroid file(s)",
            file=sys.stderr,
            flush=True,
        )
        for item, exc in failures:
            print(f"  failed: {item}: {exc}", file=sys.stderr, flush=True)
        raise SystemExit(1)

    return sorted(results, key=lambda x: x["source"])


def main():
    args = parse_args()
    args.inputs = filter_inputs(args.inputs)
    args = derive_parallelism(args)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    print_parallelism(args)

    started = time()
    results = run_pipeline(args, output_dir)

    print(f"done {len(results)} file(s) {time() - started:.1f}s")


if __name__ == "__main__":
    main()
