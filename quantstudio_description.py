# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 11:45:43 2019

@author: jmcghee
"""

import pandas 
import os

#Prior to running sample coordinates separation, the coordinate file must file formatted using this block to 
# add accessioning DTI to the coordinate file
'''
file1=pandas.read_csv('LILAC_sample_coordinates_12_2019.csv')
print(file1.columns)
file2=pandas.read_csv('Biogen_accession.csv', usecols=['Custom ID','Accesionining Number DTI'])
file2.rename(columns={'Custom ID':'Accession'}, inplace=True)
print(file2)
out=pandas.merge(file1,file2, on='Accession', how='inner')
new_col=['Accession', 'Accesionining Number DTI', 'Plate', 'Row', 'Column', 'Sample', 
             'SUBJID', 'Time Text', 'Lesional Status']
out=out[new_col]
print(out.columns)
out.to_csv('main_out.csv')'''

#Function takes in sample coordinates file and biogen file. Takes the total tield from biogen file and merges them with 
#sample coordinates file
def sample_coordinates_separattion(file1, file2):
    path='Plate output'
    #reads in biogen file to extract total yields
    data=pandas.read_csv('Biogen_accession.csv', usecols=['Custom ID','Total Yield pg'])
    #Rename columns to merge with coordinates file
    data.rename(columns={'Custom ID':'Accession'}, inplace=True)
    #Read in and merge to get the total yields
    data2=pandas.read_csv('Sample_coordinates_2_main.csv')
    df=pandas.merge(data,data2, on='Accession', how='inner')
    dff=pandas.DataFrame(df)
    columns=list(df.columns)
    
    #Rename and organize the columns
    new_col=['Accession', 'Unnamed: 0', 'Accesionining Number DTI', 'Plate', 'Row', 'Column', 'Sample', 
             'SUBJID', 'Time Text', 'Lesional Status', 'Total Yield pg']
    df=df[new_col]
    
    #extract numers from the plate number into a separate column with just numbers for sorting
    numbers=df['Plate'].str.extract('(\d+)').astype(int)
    print(numbers)
    df['sorting_nums']=numbers
    #sort by the extarcted plate numbers
    df.sort_values(by=['sorting_nums'], inplace=True)
    print(df[:20])
    
    #split the dataframe by grouping together values in the 'sorting numbers' column
    plate_id=[]
    for j, k in df.groupby(['sorting_nums']):
        plate_id.append(j)
        df=pandas.DataFrame(k)
        df.drop(['sorting_nums'], axis=1, inplace=True)
        print(df)
        #write out each 'gruped' dataframe
        df.to_csv(os.path.join(path,'Plate_'+str(j) +'.csv'))
    
    print(plate_id)
    print(columns)
    dff=dff.sort_values(by=['Total Yield pg'])
    print(dff)
    
    df.to_csv('output.csv')
 
#Takes in the ouptut files from sample_coordinates_separation and the plate array layout to output format for machine
def plate_config(f):
    path='Plate_map_out'
    f1='Plate_'+str(f)+'.csv'
    f2='slide_layouts_plate'+str(f)+'.csv'
    print(f1)
    print(f2)
    #output file from sample_coordinates_separation
    output_plate=pandas.read_csv(f1)
    #slide layout from biogen
    slide_layout=pandas.read_csv(f2, usecols=[1,2,3,4,5])
    
    #creates a dictionary so all sample names in the layout file can be converted to accesioning numbers
    dic=output_plate.set_index('Sample')['Accesionining Number DTI'].to_dict()
    print(dic)
    #Makes selected row the columns
    regen=slide_layout.replace(dic)
    regen.columns=regen.iloc[1]
    #drops unecessary rows from dataframe
    regen.drop(regen.index[[0,1]], inplace=True)
    #renames the first column as 'Lane'
    regen.rename(columns={regen.columns[0]: 'Lane'}, inplace=True)

    #transposes the regen dataframe
    regen_t=regen.T
    regen_t.reset_index(inplace=True)
    #makes selected row the column headers
    regen_t.columns=regen_t.iloc[0]
    #drops unecessary row
    regen_t.drop(regen_t.index[0], inplace=True)
    
    #creates lists for each new row for a dataframe by selecting the position values from previous dataframe
    out1=[regen_t.iloc[0,1], regen_t.iloc[0,2],regen_t.iloc[1,1],regen_t.iloc[1,2],regen_t.iloc[2,1],regen_t.iloc[2,2],
          regen_t.iloc[3,1],regen_t.iloc[3,2]]
    out2=[regen_t.iloc[0,3], regen_t.iloc[0,4],regen_t.iloc[1,3],regen_t.iloc[1,4],regen_t.iloc[2,3],regen_t.iloc[2,4],
          regen_t.iloc[3,3],regen_t.iloc[3,4]]
    out3=[regen_t.iloc[0,5], regen_t.iloc[0,6],regen_t.iloc[1,5],regen_t.iloc[1,6],regen_t.iloc[2,5],regen_t.iloc[2,6],
          regen_t.iloc[3,5],regen_t.iloc[3,6]]
    out4=[regen_t.iloc[0,7], regen_t.iloc[0,8],regen_t.iloc[1,7],regen_t.iloc[1,8],regen_t.iloc[2,7],regen_t.iloc[2,8],
          regen_t.iloc[3,7],regen_t.iloc[3,8]]
    out5=[regen_t.iloc[0,9], regen_t.iloc[0,10],regen_t.iloc[1,9],regen_t.iloc[1,10],regen_t.iloc[2,9],regen_t.iloc[2,10],
          regen_t.iloc[3,9],regen_t.iloc[3,10]]
    out6=[regen_t.iloc[0,11], regen_t.iloc[0,12],regen_t.iloc[1,11],regen_t.iloc[1,12],regen_t.iloc[2,11],regen_t.iloc[2,12],
          regen_t.iloc[3,11],regen_t.iloc[3,12]]

    #creates an empty dataframe
    output=pandas.DataFrame()
    output['Index']=['A','B','C','D','E','F','G','H']
    #adds list data into the new dataframe
    output['1']=out1; output['2']=out2; output['3']=out3; output['4']=out4; output['5']=out5; output['6']=out6
    #sets the index to the index column
    output.set_index('Index', inplace=True)
    #deletes the name/title of the index column
    del output.index.name

    #swaps the keys and values for the dictionary
    origin_dic=dict((v,k) for k,v in dic.items())
    #replaces the accession numbers with sample names
    origin=output.replace(origin_dic)
    
    #create a blank row
    missing=['']
    spacer=pandas.DataFrame(index=missing)
    #append the two dataframs with the blank datafram in between.
    appended=[output,spacer,origin]
    #concat all dataframes and name the index
    final_out=pandas.concat(appended, keys=['Accesionining Number DTI','','Sample'], sort=False)
    print(final_out)
    
    
    path='plate_map_out'
    final_out.to_csv(os.path.join(path,'output_plate'+str(f)+'.csv'))