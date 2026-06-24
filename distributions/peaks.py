# motif-quantification/distributions/peaks.py

import numpy as np


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
    array = np.asarray(array, dtype=np.float64).reshape(-1)

    if array.size == 0:
        return []

    maxes = minpoint_reduction(array, mindist)
    peakparameters = boundary_finding(maxes, array)

    return peakparameters.tolist()


def moving_average(array, width):
    array = np.asarray(array, dtype=np.float64).reshape(-1)

    if width <= 1 or array.size < width:
        return array

    kernel = np.ones(width, dtype=np.float64) / width
    return np.convolve(array, kernel, mode="same")
