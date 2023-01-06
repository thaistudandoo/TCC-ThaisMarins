from pydoc import locate
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def week_post_processing():

    #--- Generate Cortisol Plot by paper ---
    decades = ['30-40', '40-50', '50-60', '60-70', '70-80', '80-90']
    headers = ['index', 'avg','max','min','std']
    x = [30, 35, 45, 55, 65, 75, 85]

    # Cortisol
    fig, (ax1) = plt.subplots(1,1)  
    for i in range(1,6):
        fname = f'Output/cortisol/{i}_cortisol.csv'
        valor = pd.read_csv(fname, header=None)
        y = valor.T
        d = decades[i]
        label = f'Década {d}'
        x = np.linspace(0,7,7000)
        ax1.plot(x,y[0],'.',label=label)
        
    ax1.legend(bbox_to_anchor = (0.5, -0.15), loc='upper center', fontsize = 18, fancybox=True, shadow=True, ncol=5)    
    ax1.set_ylabel('Cortisol (ng/day)', fontsize = 18)
    ax1.set_xlabel('Time (days)', fontsize = 18)
    #ax1.figure(figsize=(10, 10))
    fig.set_figwidth(15) 
    fig.set_figheight(6) 
    fig.tight_layout()
    plt.savefig('Output/cortisol/cortisol.png', bbox_inches='tight')     


    # TNF
    fig, (ax2) = plt.subplots(1,1)  
    for i in range(1,6):
        fname = f'Output/tnf/{i}_TNF.csv'
        valor = pd.read_csv(fname, header=None)
        y = valor.T
        d = decades[i]
        label = f'Década {d}'
        x = np.linspace(0,7,7000)
        out_TNF = y[0]
        #out_TNF = 100 * (out_TNF - min(out_TNF)) / (max(out_TNF) - min(out_TNF))
        ax2.plot(x,out_TNF,'.',label=label)
       
        
    ax2.legend(bbox_to_anchor = (0.5, -0.15), loc='upper center', fontsize = 18, fancybox=True, shadow=True, ncol=5)    
    ax2.set_ylabel('TNF-α concentrations \n (relative values)', fontsize = 18)
    ax2.set_xlabel('Time (days)', fontsize = 18)
    #ax1.figure(figsize=(10, 10))
    fig.set_figwidth(15) 
    fig.set_figheight(6) 
    fig.tight_layout()
    plt.savefig('Output/tnf/tnf.png', bbox_inches='tight')

    # S.aureus
    fig, (ax3) = plt.subplots(1,1)  
    for i in range(1,6):
        fname = f'Output/bacteria/{i}_bacteria.csv'
        valor = pd.read_csv(fname, header=None)
        y = valor.T
        d = decades[i]
        label = f'Década {d}'
        x = np.linspace(0,7,7000)
        ax3.plot(x,y[0],'.',label=label)
        
    ax3.legend(bbox_to_anchor = (0.5, -0.15), loc='upper center', fontsize = 18, fancybox=True, shadow=True, ncol=5)    
    ax3.set_ylabel('S. aureus \n (cells/mm³)' , fontsize = 18)
    ax3.set_xlabel('Time (days)', fontsize = 18)
    fig.set_figwidth(15) 
    fig.set_figheight(6) 
    fig.tight_layout()
    plt.savefig('Output/bacteria/bacteria.png', bbox_inches='tight') 

    # Citocinas
    fig, (ax4) = plt.subplots(1,1)

    #for i in range(1,6):
    fname_il6  = f'Output/il6/6_il6.csv'
    fname_il8  = f'Output/il8/6_il8.csv'
    fname_il10 = f'Output/il10/6_il10.csv'
    valor_il6 = pd.read_csv(fname_il6, header=None)
    valor_il8 = pd.read_csv(fname_il8, header=None)
    valor_il10 = pd.read_csv(fname_il10, header=None)
    y_il6 = valor_il6.T
    y_il8 = valor_il8.T
    y_il10 = valor_il10.T
    d = decades[5]
    label_il6 =  f'IL-6  Década {d}'
    label_il8 =  f'IL-8  Década {d}'
    label_il10 = f'IL-10 Década {d}'
    x = np.linspace(0,7,7000)
    ax4.plot(x,y_il6[0],'.',label=label_il6)
    ax4.plot(x,y_il8[0],'.',label=label_il8)
    ax4.plot(x,y_il10[0],'.',label=label_il10)
        
    ax4.legend(bbox_to_anchor = (0.5, -0.15), loc='upper center', fontsize = 18, fancybox=True, shadow=True, ncol=5)    
    ax4.set_ylabel('Cytokine concentrations \n (relative values)', fontsize = 18)
    ax4.set_xlabel('Time (days)', fontsize = 18)
    fig.set_figwidth(15) 
    fig.set_figheight(6) 
    fig.tight_layout()
    plt.savefig('Output/citocinas/citocina.png', bbox_inches='tight')

   
    print('Post-processing done. Bye!')



'''
# Cytokines
     fig, (ax1) = plt.subplots(1,1)
     ax1.plot(t, out_TNF,'purple',  linewidth=3, label="TNF α")
     ax1.plot(t, out_IL6, 'b', linewidth=3,  label="IL-6")
     ax1.plot(t, out_IL8, 'r--',  linewidth=3, label="IL-8")
     ax1.plot(t, out_IL10, 'orange',  linewidth=3, label="IL-10")


     #ax1.legend( ncol = 4, bbox_to_anchor = (0.5,-0.13), loc='upper center', fontsize = 18)
     ax1.legend(bbox_to_anchor = (1,.5), loc='center left', fontsize = 18)
     ax1.set_xlabel('Time (days)', fontsize = 18)
     ax1.set_ylabel('Cytokine concentrations \n (relative values)', fontsize = 18)
     ax1.tick_params(labelsize=18)

     fig.set_figwidth(10) 
     fig.set_figheight(6) 
     fig.tight_layout()
     filename = f'{day}_Cytokines.png'
     plt.savefig(filename)

'''







'''
    for i in range(1,6):
        fname = f'Output/cortisol/1_cortisol.csv'
        valor = pd.read_csv(fname, header=None)
        y = valor.T
        d = decades[i]
        label = 'Cortisol without \nglucose influence'
        x = np.linspace(0,1,7000)
        plt.plot(x,y[0],'-',lw=3,label=label)
        plt.legend(loc='upper right',fontsize=16)
        
        plt.ylabel('Cortisol (ng/day)', fontsize = 18)
        plt.xlabel('Time (days)', fontsize = 18)
        plt.savefig('Output/cortisol/cortisol_without.png',bbox_inches='tight') 

'''

if __name__ == "__main__":
    week_post_processing()