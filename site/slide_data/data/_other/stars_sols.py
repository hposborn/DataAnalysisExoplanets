
def derive_lum(mag,dist):
    absmag = mag - 5*np.log10(dist)+5
    #assuming absmag==bolometric mag...
    lum = 10**(0.4*(4.74 - absmag)) #Finish this function!
    return lum
results['lum'] = derive_lum(results['phot_g_mean_mag'],1000/results['parallax'])
