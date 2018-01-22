''' Current functions used to calculate the scale height of the different
    Milky Way structural components. Includes reading in of current C3/C6
    input files (w/ spectro and w/o spectro parameters).
    Edited: 22/01/2018; Ben Rendle '''


import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.odr.odrpack as odrpack

def space_density(obs):
    ''' Calculate the space density for a given band above the Galactic plane '''
    # obs['Zabs'] = np.abs(obs['Z'])
    a = min(obs['Z'])
    b = max(obs['Z'])
    theta = 0.5*np.sqrt(116)
    t_theta = np.tan(theta*(np.pi/180))**2
    bins = np.linspace(np.log10(100),np.log10(5000),19)
    Z = np.log10(np.abs(obs['Z'])*1000)
    hist, bin_edges = np.histogram(Z, bins = bins)
    # print(bin_edges)
    # print(hist)
    volume = []
    # volume of square based pyramid section
    for i in range(len(bin_edges)-1):
        V1 = (4*t_theta*(10**bin_edges[i])**3)/3
        V2 = (4*t_theta*(10**bin_edges[i+1])**3)/3
        Vol = V2 - V1
        volume.append(Vol)
    # volume of cube
    # for i in range(len(bin_edges)-1):
    #     V1 = (2*np.sqrt(t_theta)*10**max(bin_edges))**2 * (10**(bin_edges[i+1]) - 10**(bin_edges[i]))
    #     volume.append(V1)
    rho = []
    rho10 = np.zeros(len(volume))
    for i in range(len(volume)):
        density = hist[i]/volume[i]
        rho10[i] = density
        rho.append(np.log10(density)) # In units number of stars/pc^3

    print(rho)
    bin_10 = 10**bin_edges[1:]
    print( len(bin_10), len(rho10))
    ''' Least squares fitting '''
    def f(Par,z):
        return Par[0]*np.exp(-z/Par[1])
    mpar, cpar, empar, ecpar = [], [], [], []
    linear = odrpack.Model(f)
    mydata = odrpack.RealData(bin_10[5:12], rho10[5:12])#, sy=df2['feh_err'])
    myodr = odrpack.ODR(mydata, linear, beta0=[1e-3, 500],maxit=20000)
    myoutput = myodr.run()
    mpar.append(myoutput.beta[0])
    cpar.append(myoutput.beta[1])
    empar.append(myoutput.sd_beta[0])
    ecpar.append(myoutput.sd_beta[1])
    # myoutput.pprint()

    b = bin_10[14:]
    b = np.delete(b,[1])
    r = rho10[14:]
    r = np.delete(r,[1])
    def f1(Par,z):
        return Par[0]*np.exp(-z/Par[1])
    mpar1, cpar1, empar1, ecpar1 = [], [], [], []
    linear1 = odrpack.Model(f1)
    mydata1 = odrpack.RealData(b, r)#, sy=df2['feh_err'])
    myodr1 = odrpack.ODR(mydata1, linear1, beta0=[5e-4, 1000],maxit=20000)
    myoutput1 = myodr1.run()
    mpar1.append(myoutput1.beta[0])
    cpar1.append(myoutput1.beta[1])
    empar1.append(myoutput1.sd_beta[0])
    ecpar1.append(myoutput1.sd_beta[1])
    # myoutput1.pprint()

    # def f2(Par,z):
    #     return Par[0]*np.exp(-z/Par[1])
    # mpar2, cpar2, empar2, ecpar2 = [], [], [], []
    # linear2 = odrpack.Model(f2)
    # mydata2 = odrpack.RealData(bin_10[4:], rho10[4:])#, sy=df2['feh_err'])
    # myodr2 = odrpack.ODR(mydata2, linear2, beta0=[1e-4, 600],maxit=20000)
    # myoutput2 = myodr2.run()
    # mpar2.append(myoutput2.beta[0])
    # cpar2.append(myoutput2.beta[1])
    # empar2.append(myoutput2.sd_beta[0])
    # ecpar2.append(myoutput2.sd_beta[1])
    # myoutput2.pprint()

    plt.figure()
    plt.scatter(bin_10,rho)
    plt.scatter(bin_10[5:12],rho[5:12],color='r')
    plt.scatter(bin_10[13:],rho[13:],color='y')
    # plt.scatter(b,np.log10(r))
    plt.plot(bin_10,np.log10(myoutput.beta[0]*np.exp(-bin_10/myoutput.beta[1])),color='orange', \
            label=r'$\rho_{0} =$ %.4g pc$^{-3}$, H$_{\rm{z}} = $ %.4g pc'%(mpar[0],cpar[0]))
    plt.plot(bin_10,np.log10(myoutput1.beta[0]*np.exp(-bin_10/myoutput1.beta[1])),color='orange',linestyle='--', \
            label=r'$\rho_{0} =$ %.4g pc$^{-3}$, H$_{\rm{z}} = $ %.4g pc'%(mpar1[0],cpar1[0]))
    plt.ylim(-7.75,-3.75)
    plt.xlim(100, 5100)
    plt.tick_params(labelsize=15)
    plt.ylabel(r'Log Space Density',fontsize=20)
    plt.xlabel(r'Z [pc]',fontsize=20)
    plt.legend(prop={'size':10})
    plt.tight_layout()
    # plt.savefig('/home/bmr135/Dropbox/K2Poles/pop_trends/211117/C3_Gilmore_Reid.png')
    # plt.show()

def space_density2(obs):
    ''' Improved calculation of the volume for determining the space density '''
    theta = np.deg2rad(np.sqrt(116))
    t1 = theta/2.0
    alpha_c3 = np.deg2rad(61.4) # degrees
    alpha_c6 = np.deg2rad(50.4) # degrees

    bins = np.linspace(np.log10(100),np.log10(5000),19)
    Z = np.log10(np.abs(obs['Z'])*1000)
    hist, bin_edges = np.histogram(Z, bins = bins)
    print(bin_edges)
    print(hist)
    volume = []

    for i in range(len(bin_edges)-1):
        x2 = lambda x: ((10**bin_edges[i])**3)*(1-np.cos(theta))*np.cos(t1)*(np.sin(alpha_c3 + t1 - x)**-3)
        x3 = lambda xa: ((10**bin_edges[i+1])**3)*(1-np.cos(theta))*np.cos(t1)*(np.sin(alpha_c3 + t1 - xa)**-3)

        vol1 = integrate.quad(x2,0,theta)
        vol2 = integrate.quad(x3,0,theta)
        Vol = vol2[0] - vol1[0]
        volume.append(Vol)

    print(volume)
    rho = []
    rho10 = np.zeros(len(volume))
    for i in range(len(volume)):
        density = hist[i]/volume[i]
        rho10[i] = density
        rho.append(np.log10(density)) # In units number of stars/pc^3

    print(rho)

    bin_10 = 10**bin_edges[1:]
    # x = pd.DataFrame()
    # x['Z'] = bin_10
    # x['log_rho'] = np.log10(rho10)
    # x['log_sig_rho'] = np.log10(np.sqrt(hist)/volume)
    # x['rho'] = 10**np.log10(rho10)
    # x['sig_rho'] = 10**np.log10(np.sqrt(hist)/volume)
    # print(x)
    # x.to_csv('/home/ben/Dropbox/K2Poles/pop_trends/Ioana_Scale_Heights/C6',index=False)
    sig_rho = np.sqrt(hist)/volume

    ''' Least squares fitting '''
    def f(Par,z):
        return Par[0]*np.exp(-z/Par[1])
    mpar, cpar, empar, ecpar = [], [], [], []
    linear = odrpack.Model(f)
    mydata = odrpack.RealData(bin_10[6:11], rho10[6:11], sy=sig_rho[6:11])
    myodr = odrpack.ODR(mydata, linear, beta0=[1e-3, 500],maxit=20000)
    myoutput = myodr.run()
    mpar.append(myoutput.beta[0])
    cpar.append(myoutput.beta[1])
    empar.append(myoutput.sd_beta[0])
    ecpar.append(myoutput.sd_beta[1])
    # myoutput.pprint()

    b = bin_10[12:]
    # b = np.delete(b,[1])
    r = rho10[12:]
    # r = np.delete(r,[1])
    sig = sig_rho[12:]
    # sig = np.delete(sig,[1])
    def f1(Par,z):
        return Par[0]*np.exp(-z/Par[1])
    mpar1, cpar1, empar1, ecpar1 = [], [], [], []
    linear1 = odrpack.Model(f1)
    mydata1 = odrpack.RealData(b, r, sy=sig)
    myodr1 = odrpack.ODR(mydata1, linear1, beta0=[5e-4, 1000],maxit=20000)
    myoutput1 = myodr1.run()
    mpar1.append(myoutput1.beta[0])
    cpar1.append(myoutput1.beta[1])
    empar1.append(myoutput1.sd_beta[0])
    ecpar1.append(myoutput1.sd_beta[1])
    # myoutput1.pprint()

    plt.figure()
    plt.scatter(bin_10,rho)
    # plt.scatter(bin_10[6:11],rho[6:11],color='r')
    # plt.scatter(bin_10[12:16],rho[12:16],color='y')
    # plt.scatter(b,np.log10(r))
    plt.plot(bin_10,np.log10(myoutput.beta[0]*np.exp(-bin_10/myoutput.beta[1])),color='orange', \
            label=r'$\rho_{0} =$ %.4g pc$^{-3}$, H$_{\rm{z}} = $ %.4g pc'%(mpar[0],cpar[0]))
    plt.plot(bin_10,np.log10(myoutput1.beta[0]*np.exp(-bin_10/myoutput1.beta[1])),color='orange',linestyle='--', \
            label=r'$\rho_{0} =$ %.4g pc$^{-3}$, H$_{\rm{z}} = $ %.4g pc'%(mpar1[0],cpar1[0]))
    plt.ylim(-7.75,-3.5)
    plt.xlim(100, 5100)
    plt.tick_params(labelsize=15)
    plt.ylabel(r'Log Space Density',fontsize=20)
    plt.xlabel(r'Z [pc]',fontsize=20)
    plt.legend(prop={'size':10})
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])


if __name__ == "__main__":

    C3_True = pd.read_csv('.../C3')
    C3_spec = pd.read_csv('.../C3_spectro')
    C6_True = pd.read_csv('.../C6')
    C6_spec = pd.read_csv('.../C6_spectro')

    space_density2(C6_True)
    space_density2(C3_True)

    plt.show()
