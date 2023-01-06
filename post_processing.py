from pydoc import locate
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def post_processing():

    #--- Generate Cortisol Plot by paper ---
    decades = ['30-40', '40-50', '50-60', '60-70', '70-80', '80-90']
    headers = ['index', 'avg','max','min','std']
    x = [30, 35, 45, 55, 65, 75, 85]

    # (Female) Read csv file to plot
    valor_f = pd.read_csv('average_cortisol_female.csv', names = headers)
    #print(valor_f)
    y_f = valor_f['avg']
    #print("tamanho yf:", np.size(y_f))
    #print("y_f:", y_f)
    std_f = valor_f['std']
    max_f = valor_f['max']
    min_f = valor_f['min']

    # (Male) Read csv file to plot
    valor_m = pd.read_csv('average_cortisol_male.csv', names = headers)
    #print(valor_m)
    y_m = valor_m['avg']
    std_m = valor_m['std']
    max_m = valor_m['max']
    min_m = valor_m['min']
  
    headers = ['decade', 'value']
    cortisol_female = pd.read_csv('cortisol_data_female.csv', names = headers)
    cortisol_female_5 = cortisol_female.tail(6)
    cortisol_male = pd.read_csv('cortisol_data_male.csv', names = headers)
    cortisol_male_5 = cortisol_male.tail(6)
    
    plt.figure()
    ##plt.plot(x,y, label="Cortisol", c= 'g', linewidth=0.75)
    #plt.plot(x,y, label="Cortisol", c= 'g', linewidth=2)
    plt.errorbar(x, y_f, std_f, fmt='ok', lw=1, label='Average (Model)')
    #plt.errorbar(y_f['index'], y_f, std_f, fmt='ok', lw=1, label='Average (Model)')
    plt.errorbar(x,y_f,[y_f-min_f, max_f-y_f],fmt='.k', ecolor='gray', lw=1)
    plt.plot(cortisol_female_5['decade'],cortisol_female_5['value'],'or',lw=5, label='Exp. Data')
    plt.legend (fontsize = 18,loc='upper center',bbox_to_anchor=(0.5, 1.2), fancybox=True, shadow=True)
    plt.xlabel('Age (years)', fontsize = 18)
    plt.ylabel('Cortisol (ng/day)', fontsize = 18)
    plt.savefig('Cortisol_box_female.png',bbox_inches='tight')
    
    plt.figure()
    ##plt.plot(x,y, label="Cortisol", c= 'g', linewidth=0.75)
    #plt.plot(x,y, label="Cortisol", c= 'g', linewidth=2)
    plt.errorbar(x,y_m,std_m,fmt='sk', lw=1,label='Average (Model)')
    plt.errorbar(x,y_m,[y_m-min_m, max_m-y_m],fmt='.k', ecolor='gray', lw=1)
    plt.plot(cortisol_male_5['decade'],cortisol_male_5['value'],'ob',lw=5, label='Exp. Data')
    plt.legend (fontsize = 18,loc='upper center',bbox_to_anchor=(0.5, 1.2), fancybox=True, shadow=True)
    plt.xlabel('Age (years)', fontsize = 18)
    plt.ylabel('Cortisol (ng/day)', fontsize = 18)
    plt.savefig('Cortisol_box_male.png',bbox_inches='tight')
    
    plt.figure()
    #plt.plot(x,max_f,'r',lw=2,label='Model Female')
    #plt.plot(x,max_m,'b',lw=2, label='Model Male')
    plt.plot(cortisol_female_5['decade'],cortisol_female_5['value'],'-r',markersize=10, label='Model Female')
    plt.plot(cortisol_male_5['decade'],cortisol_male_5['value'],'-b',markersize=10, label='Model Male')
    plt.plot(cortisol_female_5['decade'],cortisol_female_5['value'],'or',markersize=10, label='Exp. Data Female')
    plt.plot(cortisol_male_5['decade'],cortisol_male_5['value'],'sb',markersize=10, label='Exp. Data Male')
    plt.xlabel('Age (years)', fontsize = 18)
    plt.legend(fontsize=16)
    plt.ylabel('Cortisol (ng/day)', fontsize = 18)
    plt.savefig('Cortisol_max.png',bbox_inches='tight')
    
   
    ### activated macrophages #################################
    plt.figure()
    for i in range(0,5):
        fname = f'Output/female/{i}_ma.csv'
        valor = pd.read_csv(fname, header=None)
        y = valor.T
        d = decades[i]
        label = f'Decade {d}'
        x = np.linspace(0,1,1000)
        plt.plot(x,y[0],'.',label=label)
        plt.legend()

    plt.ylabel('Activated Macrophages concentrations', fontsize = 18)
    plt.xlabel('Time (days)', fontsize = 18)
    plt.savefig('ma.png',bbox_inches='tight')

    ### TNF
    plt.figure()
    for i in range(1,6):
        fname = f'Output/female/{i}_TNF.csv'
        valor = pd.read_csv(fname, header=None)
        y = valor.T
        d = decades[i]
        label = f'Decade {d}'
        x = np.linspace(0,1,1000)
        out_TNF = y[0]
        out_TNF = 100 * (out_TNF - min(out_TNF)) / (max(out_TNF) - min(out_TNF))
        plt.plot(x,out_TNF,'.',label=label)
        plt.legend(fontsize=16)
        
    plt.ylabel('TNF-α concentrations \n (relative values)', fontsize = 18)
    plt.xlabel('Time (days)', fontsize = 18)
    plt.savefig('tnf-alpha.png',bbox_inches='tight') 
    
    ### Cortisol without glucose
    
    plt.figure()
    fname = f'Output/female/0_cortisol.csv'
    valor = pd.read_csv(fname, header=None)
    y = valor.T
    d = decades[i]
    label = 'Cortisol without \nglucose influence'
    x = np.linspace(0,1,1000)
    plt.plot(x,y[0],'-',lw=3,label=label)
    plt.legend(loc='upper right',fontsize=16)
    
    plt.ylabel('Cortisol (ng/day)', fontsize = 18)
    plt.xlabel('Time (days)', fontsize = 18)
    plt.savefig('cortisol_without.png',bbox_inches='tight') 
    
   
    
    ### Cortisol
    plt.figure()
    for i in range(1,6):
        fname = f'Output/female/{i}_cortisol.csv'
        valor = pd.read_csv(fname, header=None)
        y = valor.T
        d = decades[i]
        label = f'Decade {d}'
        x = np.linspace(0,1,1000)
        plt.plot(x,y[0],'.',label=label)
        
    plt.legend(fontsize=16)    
    plt.ylabel('Cortisol (ng/day)', fontsize = 18)
    plt.xlabel('Time (days)', fontsize = 18)
    plt.savefig('cortisol.png',bbox_inches='tight') 
    
    # Glucose
    plt.figure()
    fname = f'Output/female/0_glucose.csv'
    valor = pd.read_csv(fname, header=None)
    y = valor.T
    d = decades[i]
    label = 'Glucose intake'
    x = np.linspace(0,12,720000)
    plt.plot(x,y[0],'-',lw=3,label=label)
    
    plt.legend(loc='upper right',fontsize=16)
    plt.ylabel('Glucose (mmol)', fontsize = 18)
    plt.xlabel('Time (hours)', fontsize = 18)
    plt.savefig('glucose.png',bbox_inches='tight') 
    
    
    print('Post-processing done. Bye!')

 
if __name__ == "__main__":
    post_processing()