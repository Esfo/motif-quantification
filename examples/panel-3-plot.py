for status, distgroups in plotdicts.items():
    for distkey, dists in distgroups.items():
        chargeorder = sorted(dists)
        chlen = len(dists)
        cf = pd.DataFrame()
        chargefigures = {c:n for n, c in enumerate(chargeorder)}
        fig, ax = plt.subplots(ncols=chlen, nrows=8, figsize=(6,8), sharex='col', sharey='row')
        fig.subplots_adjust(hspace=0.05, wspace=0.05)
        cg = list(dists.values())
        chargedistbounds = [np.inf, 0]
        intensitylineup = [distributionintensities[g] for g in cg]
        flatintensities = itertools.chain(*intensitylineup)
        maxmain = max(flatintensities)
        masslineup = [distributionmasses[g]*distributioncharges[g]-proton*distributioncharges[g] for g in cg]
        masslineup = sorted(masslineup, key=lambda x: -x.size)
        intensitylineup = sorted(intensitylineup, key=lambda x: -x.size)
        arraysizes = [i.size for i in masslineup]
        if len(set(arraysizes)) == 1:
            arraymeans = np.array(masslineup).mean(axis=0)
            intensitysums = np.array(intensitylineup).sum(axis=0)
        else:
            matrixmax = max(arraysizes)
            arraysums = np.zeros(matrixmax)
            arraydividends = np.zeros(matrixmax)
            intensitysums = np.zeros(matrixmax)
            for n, (a, i) in enumerate(zip(masslineup, intensitylineup)):
                orderind = np.abs(masslineup[0] - a[0]).argmin()
                asize = a.size
                if asize == matrixmax:
                    arraysums += a
                    arraydividends += 1
                    intensitysums += i
                else:
                    arraysums[orderind:orderind+asize] += a
                    arraydividends[orderind:orderind+asize] += 1
                    intensitysums[orderind:orderind+asize] += i
            arraymeans = arraysums / arraydividends
        tbarwidth = 1 / len(cg)
        tspace = 0.5
        cols = dp.get_colors(len(cg))
        cmin = np.inf
        cmax = 0
        cratiolist = []
        ccbounds = [np.inf, 0]
        for n, (charge, g) in enumerate(dists.items()):
            lines = linesofdistributions[g]
            cintensities = distributionintensities[g]
            cratios = cintensities[:-1] / cintensities[1:]
            cratios[cratios < 1] = -1 / cratios[cratios < 1]
            cratiolist.append(cratios)
            abcratios = [abs(i) for i in cratios]
            cmasses = distributionmasses[g]
            if cintensities.min() < cmin:
                cmin = cintensities.min()
            if cintensities.max() > cmax:
                cmax = cintensities.max()
            concharge = distributioncharges[g]
            expdiff = proton / concharge
            basemasses = cmasses * concharge - proton * concharge
            if basemasses.size > arraymeans.size:
                orderind = np.abs(basemasses - arraymeans[0]).argmin()
                matrixmin = arraymeans.size
                meanbasediff = arraymeans - basemasses[orderind:orderind+matrixmin]
                meanbaseppms = (meanbasediff / basemasses[orderind:orderind+matrixmin]) * 1000000
                cmx = cmasses[orderind:orderind+matrixmin]
            else:
                orderind = np.abs(arraymeans - basemasses[0]).argmin()
                matrixmin = basemasses.size
                meanbasediff = arraymeans[orderind:orderind+matrixmin] - basemasses
                meanbaseppms = (meanbasediff / basemasses) * 1000000
                cmx = cmasses
            acdiffs = expdiff - np.diff(cmasses)
            basediffs = acdiffs * concharge
            conhax = chargefigures[concharge]
            cwidth = 0.5 *  len(lines)
            chargelengthextra = proton / charge * 2
            nst = regions[lines,2].min() - 0.5
            net = regions[lines,3].max() + 0.5
            nlmb = regions[lines,0].min() - chargelengthextra
            numb = regions[lines,1].max() + chargelengthextra
            nboundrec = [nlmb, numb, nst, net]
            nplotkeys = arg_coord_rectangle_overlap(nboundrec, regions[:,:4]).tolist()
            for p in nplotkeys:
                if p not in lines:
                    a = trackedgroups[p]
                    creg = regions[p]
                    ax[0][conhax].plot(a[0], a[1], '.', color='white', markersize=0.8, alpha=0.3)
                    ax[0][conhax].plot(a[0], a[1], '-', color='white', linewidth=0.4, alpha=1)
                    ax[1][conhax].plot([creg[7], creg[7]], [0, creg[5]], '-', color='white', alpha=0.5, linewidth=cwidth)
            ax[5][conhax].hlines(0, cmasses.min(), cmasses.max(), color='black', linewidth=0.3)
            ax[7][conhax].hlines(0, cmasses.min(), cmasses.max(), color='black', linewidth=0.3)
            ax[2][conhax].hlines(proton, cmasses.min(), cmasses.max(), color='black', linewidth=0.3)
            for cline in lines:
                creg = regions[cline]
                ax[1][conhax].plot([creg[7], creg[7]], [0, creg[5]], '-', color=cols[n], alpha=1, linewidth=cwidth)
                a = trackedgroups[cline]
                ax[0][conhax].plot(a[0], a[1], '.', color=cols[n], markersize=0.8, alpha=0.3)
                ax[0][conhax].plot(a[0], a[1], '-', color=cols[n], linewidth=0.4, alpha=1)
                #ax[2][conhax].plot([creg[7], creg[7]], [0, creg[4]], '-', color=cols[n], alpha=1, linewidth=cwidth)
            basenorms = basediffs / proton
            abasediffs = np.abs(basenorms / proton - basediffs)
            diffgen = 0.05
            bn = 0
            bw = 0.02
            adjacentx = cmasses[:-1] + np.diff(cmasses) / 2
            ax[5][conhax].bar(adjacentx, cratios, width=diffgen, color=cols[n], alpha=1)
            #ax[6][conhax].bar(adjacentx, cpointratios, width=diffgen, color=cols[n], alpha=1)
            chargedists = np.diff(cmasses) * charge
            ax[2][conhax].bar(adjacentx, chargedists, width=diffgen, color=cols[n], alpha=1)
            intensitypercs = cintensities[:intensitysums.size] / intensitysums[:cintensities.size]
            for ip in chargedists:
                if ip < chargedistbounds[0]:
                    chargedistbounds[0] = ip
                if ip > chargedistbounds[1]:
                    chargedistbounds[1] = ip
            ax[4][conhax].bar(cmasses, intensitypercs, width=diffgen, color=cols[n], alpha=1)
            ax[6][conhax].bar(cmx, meanbaseppms, width=diffgen, color=cols[n], alpha=1)
            for nc, cn in enumerate(cg):
                cnmasses = distributionmasses[cn]
                cncharge = distributioncharges[cn]
                cnbases = cnmasses * cncharge - proton * cncharge
                cnintensities = distributionintensities[cn]
                if basemasses.size > cnbases.size:
                    orderind = np.abs(basemasses - cnbases[0]).argmin()
                    matrixmin = cnbases.size
                    cnratio = cnintensities / cintensities[orderind:orderind+matrixmin]
                    bx = cmasses[orderind:orderind+matrixmin]+(diffgen*nc)
                else:
                    orderind = np.abs(cnbases - basemasses[0]).argmin()
                    matrixmin = basemasses.size
                    cnratio = cnintensities[orderind:orderind+matrixmin] / cintensities
                    bx = cmasses+(diffgen*nc)
                cnratiobar = np.abs(cnratio.mean() - cnratio).mean()
                if cnratiobar > ccbounds[1]:
                    ccbounds[1] = cnratiobar
                if cnratiobar < ccbounds[0] and cnratiobar > 0:
                    ccbounds[0] = cnratiobar
                cx = chargefigures[distributioncharges[cn]]
                cf.loc[g, cn] = cnratiobar
                ax[3][conhax].bar(bx, cnratio, width=diffgen, color=cols[nc], alpha=1)
                #ax[4][conhax].bar(bx, pointratio, width=diffgen, color=cols[nc], alpha=1)
            for nc, cn in enumerate(cg):
                if cn != cg:
                    maincharge = distributioncharges[cn]
                    mainmasses = distributionmasses[cn]
                    mainbasemasses = mainmasses * maincharge - proton * maincharge
                    if basemasses.size > mainbasemasses.size:
                        orderind = np.abs(basemasses - mainbasemasses[0]).argmin()
                        matrixmin = mainbasemasses.size
                        maindiffs = mainbasemasses - basemasses[orderind:orderind+matrixmin]
                        mainppm = (maindiffs / mainbasemasses) * 1000000
                        bx = cmasses[orderind:orderind+matrixmin]+(diffgen*nc)
                    else:
                        orderind = np.abs(mainbasemasses - basemasses[0]).argmin()
                        matrixmin = basemasses.size
                        maindiffs = mainbasemasses[orderind:orderind+matrixmin] - basemasses
                        mainppm = (maindiffs / mainbasemasses[orderind:orderind+matrixmin]) * 1000000
                        bx = cmasses+(diffgen*nc)
                    diffbar = np.abs(mainppm.mean() - mainppm).mean()
                    cx = chargefigures[distributioncharges[cn]]
                    #ax[5][conhax].bar(cmasses+bn, maindiffs, width=bw, color=cols[nc], alpha=1)
                    ax[7][conhax].bar(bx+bn, mainppm, width=bw, color=cols[nc], alpha=1)
                    bn += bw + bw / 2
            ax[0][conhax].set_title(''.join((str(concharge), '(', str(g), ')')), fontsize=12)
        ax[1][0].set_yscale('log')
        #ax[2][0].set_yscale('log')
        ax[2][0].set_ylim(chargedistbounds[0]*0.95, chargedistbounds[1]*1.05)
        ax[3][0].set_yscale('log')
        ax[5][0].set_yscale('symlog')
        ax[4][0].set_yscale('log')
        #ax[6][0].set_yscale('symlog')
        ax[1][0].set_ylim(cmin/2, cmax)
        ax[1][0].set_ylabel('peak area')
        ax[0][0].set_ylabel('retention time', fontsize=6)
        ax[3][0].set_ylabel('cross-charge', fontsize=7)
        ax[5][0].set_ylabel('adjacency')
        ax[7][0].set_ylabel('ppm error')
        ax[2][0].set_ylabel('charge distances', fontsize=6)
        ax[4][0].set_ylabel('intensity sum %', fontsize=6)
        ax[6][0].set_ylabel('ppm to mean', fontsize=7)
        for ch, hax in chargefigures.items():
            ax[-1][hax].tick_params(axis='x', labelrotation=-45)
            if hax == 0:
                #invisible right splines
                ax[0][hax].spines.right.set_visible(False)
                ax[1][hax].spines.right.set_visible(False)
            elif hax == chlen-1:
                #invisible left splines
                ax[0][hax].spines.left.set_visible(False)
                ax[1][hax].spines.left.set_visible(False)
                for tick in ax[0][hax].yaxis.get_major_ticks():
                    tick.tick1line.set_visible(False)
                    tick.tick2line.set_visible(False)
                for tick in ax[1][hax].yaxis.get_majorticklines():
                    tick.set_visible(False)
                for tick in ax[1][hax].yaxis.get_minorticklines():
                    tick.set_visible(False)
            else:
                #left and right invisible
                ax[0][hax].spines.left.set_visible(False)
                ax[0][hax].spines.right.set_visible(False)
                ax[1][hax].spines.left.set_visible(False)
                ax[1][hax].spines.right.set_visible(False)
                for tick in ax[0][hax].yaxis.get_major_ticks():
                    tick.tick1line.set_visible(False)
                    tick.tick2line.set_visible(False)
                for tick in ax[1][hax].yaxis.get_majorticklines():
                    tick.set_visible(False)
                for tick in ax[1][hax].yaxis.get_minorticklines():
                    tick.set_visible(False)
        supt = ': '.join((str(distkey), status))
        plt.suptitle(supt, y=0.92)
        plt.show()
        fig.clf()
        plt.close()
        gc.collect()
