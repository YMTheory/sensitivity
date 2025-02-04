import numpy as np

def gauss(x, mu, sigma):
    return 1 / np.sqrt(2*np.pi*sigma**2) * np.exp(-(x-mu)**2/2/sigma**2)

def convolve_energy_resolution(edges, conts, A):
    # The energy resolution formula used here is sigma/E = A/sqrt(E)
    # The original bin edges and conts ( N_edges = N_conts + 1)
    # This is only for uniform binning for now...
    
    # The input edges+conts combination is the fine binning
    # 1 for fine binning, 2 for coarse binning, bin width a factor time
    factor = 10
    width1 = edges[1] - edges[0]
    width2 = width1 * factor
    conts_smeared = np.zeros(len(conts))
    n_bin1 = len(edges) - 1
    n_bin2 = int(n_bin1 / factor) 
    xlow, xhig = edges[0], edges[-1]

    #print(f'The original {n_bin1} bins are from {xlow} to {xhig} with {width1} width. \n')
    
    # loop in the original fine binning array
    for i in range(n_bin1):
        centi = edges[i] + width1/2.
        conti = conts[i]
        sigmai = A * np.sqrt(centi)
        xlowi, xhigi = centi - 5*sigmai, centi + 5*sigmai
        xlowid, xhigid = int((xlowi-xlow)/width1), int((xhigi-xlow)/width1)
        tmp, prob = 0, 0
        for j in range(xlowid, xhigid):
            #centj = edges[j] + width1/2.
            centj = xlow + width1 * j + width1/2.
            tmp += conti * gauss(centj, centi, sigmai) * width1
            prob += gauss(centj, centi, sigmai) * width1
            #print(i, centi, j, centj)
            if j < 0:
                j = 0
            elif j >= n_bin1:
                j = n_bin1 - 1 
            conts_smeared[j] += conti * gauss(centj, centi, sigmai) * width1
        #if prob > 1:
        #    print(f"{i}-th bin running from {xlowid} to {xhigid}")
        #    print(f"{i}-th bin: {conti:.4f} -> {tmp:.4e}, {prob:.4f} with {conti < tmp}")
            
    
    ## The new smeared fine binning is the combination of edges + conts_smeared.
    ## Coarse binning
    edges_coarse = []
    for i in range(n_bin2):
        edges_coarse.append( xlow + width2 * i )
    conts_coarse = np.zeros(n_bin2)
    edges_coarse.append( edges_coarse[-1] + width2)
    edges_coarse = np.array(edges_coarse)
    for j in range(n_bin1):
        k = int(j/factor)
        conts_coarse[k] += conts_smeared[j]
    
    return edges_coarse, conts_coarse
        

        



