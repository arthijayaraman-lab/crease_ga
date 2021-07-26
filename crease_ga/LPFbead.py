#calculates spherical form factor
import numpy as np    

def LPFbead(qrange, sigmabead):
    import numpy as np    
    R=np.true_divide(sigmabead,2)
    QR=np.multiply(qrange,R)
    Fqb=np.multiply(np.true_divide(np.sin(QR)-np.multiply(QR,np.cos(QR)),np.power(QR,3)),3)  

    return Fqb

equivalent='''

    R=np.true_divide(sigmabead,2)
    QR=np.multiply(qrange,R)
    sinqr_qr2=np.true_divide(np.sin(QR),np.power(QR,2))
    cos_qr=np.true_divide(np.cos(QR),QR)
    Fqb=(np.multiply(3,np.true_divide(np.subtract(sinqr_qr2,cos_qr),QR))) #1=3\
    
    R=np.true_divide(sigmabead,2)
    QR=np.multiply(qrange,R)
    Fqb=np.multiply(np.true_divide(np.sin(QR)-np.multiply(QR,np.cos(QR)),np.power(QR,3)),3)
'''


Random_tests='''
    print Fqb[-1]
    Fqb=np.multiply(Fqb,Fqb)
    figsize=(4,4)
    fig, ax = plt.subplots(figsize=(figsize))
    ax.plot(qrange,Fqb,color='k',linestyle='-',ms=8,linewidth=1.3,marker='o')
    ax.plot((qrange[0],qrange[-1]),(0,0),color='k',linestyle='-',ms=8,linewidth=1.3,marker='o')
    #plt.xlim(0.003,0.2)
    #plt.ylim(2*10**(-5),20)
    plt.xlabel(r'q$\sigma$',fontsize=20)
    plt.ylabel(r'$I$(q)',fontsize=20)
    #plt.title('Chainsize distribution eps='+str(graphepsdist),fontsize=20)
    #ax.set_xscale("log")#, nonposy='clip')
    #ax.set_yscale("log")
    #ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    #ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    #ax.set_ylim(10**(-4),1.1)
    #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    #plt.tick_params(axis='both', which='both',labelsize=20)
    #plt.savefig('test_Itot.png', bbox_inches='tight')
    plt.show()
    
    
    figsize=(4,4)
    fig, ax = plt.subplots(figsize=(figsize))
    ax.plot(qrange,Fqb,color='k',linestyle='-',ms=8,linewidth=1.3,marker='o')
    ax.plot((qrange[0],qrange[-1]),(0,0),color='k',linestyle='-',ms=8,linewidth=1.3,marker='o')
    #plt.xlim(0.003,0.2)
    #plt.ylim(2*10**(-5),20)
    plt.xlabel(r'q$\sigma$',fontsize=20)
    plt.ylabel(r'$I$(q)',fontsize=20)
    #plt.title('Chainsize distribution eps='+str(graphepsdist),fontsize=20)
    ax.set_xscale("log")#, nonposy='clip')
    ax.set_yscale("log")
    #ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    #ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    #ax.set_ylim(10**(-4),1.1)
    #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    #plt.tick_params(axis='both', which='both',labelsize=20)
    #plt.savefig('test_Itot.png', bbox_inches='tight')
    plt.show()
    
    print asdfasdfasdf
    
'''