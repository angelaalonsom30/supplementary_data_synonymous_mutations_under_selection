#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import itertools

#readcodon and aa (leer_codon)
def leer_codon(codon,def_aminoacidos):
    #print(codon)
    return def_aminoacidos[def_aminoacidos[0]==codon].iloc[0][1]

#calc_codon_n, Function that uses codon number (n), the funcion gets as input the codon number(n), the csv dataframe with all substitution profiles (the 12 possible base substitutions raw counts) 
#at each of the three codon positions,  
#and the dictionary for each codon. 

def calc_codon_n(n,df_prob,dict_llenar):
    #list(dict_llenar.keys())
    for j in list(dict_llenar.keys()):
        #print(j,end=', ')
        prod = 1
        pos_en_codon = 3*n-2
        for i in j:
            prod = prod * df_prob.loc[pos_en_codon][i]
            pos_en_codon +=1
        dict_llenar[j] += [prod]
    return dict_llenar

#Main function converts raw nucleotide counts into normalized probabilities per codon position accounting for order in the sequence
def ProbsCodonAmino(file_name):
    #Reads the codon-amino-acid table and the raw substitution-counts csv
    
    print("lectura de la tabla codones_y_aminoacidos.txt")
    def_aminoacidos = pd.read_csv('codones_y_aminoacidos.txt', sep='\t',header = None)
    
    
    print('lectura de la tabla:', file_name+'.csv')
    df = pd.read_csv(file_name+'.csv')
    
    #Normalizing to transform counts to probabilities
    df.set_index('position', inplace=True)
    df['suma'] = df.sum(axis=1)
    for i in df.columns:
        df[i+'_norm'] = df[i]/df['suma']
        df.drop(i,axis=1,inplace=True)
    df.drop(['suma_norm'],axis=1,inplace=True)
    
    #renaming columns 
    df.columns = [i[0:3] for i in df.columns]
    
    #creating of a list with all ancestral-substitution codon combinations
    lista_anc_cod = list(df.columns)
    comb_anc_cod = list(itertools.product(lista_anc_cod,repeat=3))
    #len(comb_anc_cod)
    
    #Compute probabilities per codon position
    data_dict = {k: [] for k in comb_anc_cod}
    
    total_codones = df.shape[0]//3
    indice_codones = []
    print('Calculando probabilidades por cada codon:')
    for i in range(df.shape[0]//3):
        print("codon",i+1,end="\r")
        data_dict = calc_codon_n(i+1,df,data_dict)
        indice_codones.append(i+1)
        
    #Relabel keys with the resulting amino acid
    for i in comb_anc_cod:
        str_col = str(i).replace('(','').replace(')','').replace("'",'')#.replace(' ','').replace(',','')
        aminoacido = leer_codon(str_col[2] + str_col[7] + str_col[12],def_aminoacidos)
        data_dict[str_col + ' - ' + aminoacido] = data_dict.pop(i)
    
    df_codon = pd.DataFrame(data_dict,index=indice_codones)
    
    #Build and save the detailed table with all ancestral substitution codn(triplet) probabilities .csv
    df_codon.to_csv(file_name+'_codones_ancestros.csv')
    print('se guardo la tabla de probabilidades con ancestros en el archivo:',file_name+'_codones_ancestros.csv')
    
    
    comb_letras = list(itertools.product(['T','C','A','G',],repeat=3))
    
    #Marginalize out the ancestral information
    dict_3letras = {}
    for j in comb_letras:
        lista_probs = []
        cont=1
        for i in comb_anc_cod:
            if ((i[0][-1] == j[0]) & (i[1][-1] == j[1]) & (i[2][-1] == j[2])):
                str_col = str(i).replace('(','').replace(')','').replace("'",'')#.replace(' ','')#.replace(',','|')
                aminoacido = leer_codon(str_col[2] + str_col[7] + str_col[12],def_aminoacidos)
                lista_probs.append(str_col + ' - ' + aminoacido)
        str_col = str(j).replace('(','').replace(')','').replace("'",'').replace(' ','').replace(',','')
        aminoacido = leer_codon(str_col,def_aminoacidos);
        dict_3letras[str_col + ' - ' + aminoacido] = list(df_codon[lista_probs].sum(axis=1))
    
    #Build and save the collapsed table derived by summing over the ancestral-substitution-aware table
    df_3letras = pd.DataFrame(dict_3letras,index=indice_codones)

    df_3letras.to_csv(file_name+'_codones_sin_ancestros.csv')
    print('se guardo la tabla de probabilidades sin ancestros en el archivo:',file_name+'_codones_sin_ancestros.csv')

    
if __name__ == "__main__": 
    ProbsCodonAmino()