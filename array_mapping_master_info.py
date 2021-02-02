# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 15:37:02 2019

@author: jmcghee
"""

import pandas
import numpy as np
import os
import math
import matplotlib.pyplot as plt


# Pre-processes the raw input files from QuantStudio. Takes in raw data files and gene_list data
# Outputs dataframe will triplicates analyzed and worst dropped from the original raw file
def pre_processing(file1, file2):
     # read in gene list to use the gene symbols to merge with sample files
    gene_data=pandas.read_csv(file2, usecols=[0,1,2])
    # Replaces the column name in gene_list df to merge with raw data file
    gene_data.rename(columns={'Taqman_ID':'Target Name'}, inplace=True)
    # strips extension off raw input file
    sample_name=file1.replace('.txt','')
    # read in raw sample text file and skip excess information
    data=pandas.read_csv(file1, delimiter='\t', skiprows=19)[:-4]
    
    #Create new dataframe merging the gene symbols and 'type' (housekeeping vs. non-housekeeping) with sample file
    df=pandas.merge(data, gene_data, on='Target Name', how='outer')
    # Re-arange columns
    columns=['Well', 'Well Position', 'Omit', 'Sample Name','Gene Symbol','Type', 'Target Name', 'Task', 'Reporter', 'Quencher', 
           'Ct', 'Ct Mean', 'Ct SD', 'Amp Score', 'Cq Conf', 'CRTAMPLITUDE', 'HIGHSD', 'CRTNOISE', 'ROX Signal']
    df=df[columns]
    dff=pandas.DataFrame(df)
    
    # create dictionary to replace input control names with sample control names
    '''
    dic={'Control RNA rep  1|Control RNA rep  1':'HEK control','Control RNA rep  2|Control RNA rep  2':'HEK control',
         'Control RNA rep  3|Control RNA rep  3':'HEK control','Control RNA rep  4|Control RNA rep  4':'UHR control',
         'Control RNA rep  5|Control RNA rep  5':'UHR control','Control RNA rep  6|Control RNA rep  6':'UHR control'}
    '''

    # sort by gene suymbol and then by sample name
    dff=dff.sort_values(by=['Gene Symbol','Sample Name'])
    # repalce any 'Undetermined values with NaN
    dff=dff.replace(['Undetermined'],[np.nan])
    #replace all na values in the 'type' column with 'non-housekeeping
    dff['Type'].replace(np.nan,'non housekeeping', inplace=True)
    #change the data type to numeric for mean function
    dff['Ct']=pandas.to_numeric(dff['Ct'])
    print(dff)

    #replace all control input names with values from the dictionary
    '''
    dff['Sample Name'].replace(dic, inplace=True)
    '''
    # create a drop list that will contain index numbers to be dropped
    drop_list=[]
    # loop through grouped dataframes
    for v, r in dff.groupby(['Sample Name','Gene Symbol']):
        count_list=[]
        # Length setting will skip duplicates with more than 3
        if len(r)==3:
            # calculate mean Ct value for each sliced dataframe
            mean=r['Ct'].mean()
            # loop through each Ct value in the sliced dataframe
            for i in r['Ct']:
                # Get the absolute difference between Ct values for each sample and the df mean
                x=abs(i-mean)
                # append each absolute difference to a list
                count_list.append(x)
        
        # Repalce any nan values in the list 
        train = [10 if math.isnan(x) else x for x in count_list]
        # get the largest number in the list, will be the outlier or a nan value
        largest=max(train)
        # will get the index for the largest number in the sliced df list
        number=train.index(largest)
        # will append the index number from the dataframe using the index number from the sliced list
        drop_list.append(r.index[number])
        
    # create a new dataframe dropping all the indexed samples from triplicates
    new_df=dff.drop(drop_list)
    
    std=[]
    # loop through the dataframe, now duplicates, by sliced df
    for k, j in new_df.groupby(['Sample Name','Gene Symbol']):
        # calculate the std deviation for each sliced df
        standards=j["Ct"].std()
        # Take an index number for each sliced df and the std deviation into a list
        list_a=[j.index[0], standards]
        # append this list to the main list creating a list of lists
        std.append(list_a)
    
    #create a list of headers
    columns_headers=['Index','New Ct SD']
    # create a new dataframe for the std deviations
    out_df=pandas.DataFrame(data=std, columns=columns_headers)
    # Set the index for the std deviation list to index numbers
    out_df.set_index(['Index'], inplace=True)
    # Create the output df joining the 'duplicate df' and the std deviation df
    pre_output=new_df.join(out_df,how='outer')
    
    # reorganize the columns
    out_columns=['Well', 'Well Position', 'Omit', 'Sample Name','Gene Symbol','Type', 'Target Name', 'Task', 'Reporter', 'Quencher', 
           'Ct', 'Ct Mean', 'Ct SD','New Ct SD', 'Amp Score', 'Cq Conf', 'CRTAMPLITUDE', 'HIGHSD', 'CRTNOISE', 'ROX Signal']
    pre_output=pre_output[out_columns]
    pre_output.to_csv(sample_name+'_processed.txt', sep='\t', index=False)
    # Outputs processed raw data file with triplicates analyzed and updated std deviation calculated
    print(pre_output)

    
#pre_processing('_export.txt', 'gene_list.csv')



#Takes in raw data file from QuantStudio and filters out bad samples
def array_filtering(file):
    #reads in raw quantstudio file
    file_name=file.replace('.txt','')
    data=pandas.read_csv(file, delimiter='\t', skiprows=19)[:-4]
    cols=data.columns
    #creates path for bad_samples output
    path='filtered_samples'
    print(cols)
    
    #loops through the input file and extracts bad samples to new list and then drops them from the dataframe
    bad_samples=[]
    for q, r in data.iterrows():
        #Check for controls and either pass or continue to filter out control samples
        if 'Control' in r['Sample Name']:
            pass
        #set parameters
        if r['Amp Score'] < 1.23 or r['Ct SD'] >0.50:
            bad_samples.append(r)
            #drops bad samples from dataframe after putting them into separate list
            data.drop(q, inplace=True)
    #takes bad sample list and puts them into a dataframe
    bad_df=pandas.DataFrame(bad_samples)
    bad_df.to_csv(os.path.join(path,file_name+'_bad_samples.csv'))
    data.to_csv(file_name+'_filtered.txt', sep='\t', index=False )




#Takes in raw sample file after processing and merges duplicates, contains delta Ct values, and housekeeping values
def array_mapping(pre_output):
    #read in file
    data=pandas.read_csv(pre_output, delimiter='\t')
    print(data)
    # strips extension off raw input file
    sample_name=pre_output.replace('.txt','')
    #group all triplicates togther and calculate the mean for Ct values
    samples=data.groupby(['Sample Name', 'Gene Symbol','Target Name','Type'],as_index=False, sort=False).mean()
    samples.rename(columns={'Ct':'New Ct mean'}, inplace=True)
    #separate samples into housekeeping and non-housekeeping dataframes
    house_genes=samples[samples['Type']=='housekeeping']
    array_samples=samples[samples['Type']=='non housekeeping']
    
    perc_cv=[]
    #looping through the housekeeping dataframe and calculating the percent cv
    for f, row in house_genes.iterrows():
        std_dev=row['Ct SD']; cv_mean=row['Ct Mean']
        percent_cv=(std_dev/cv_mean)*100
        perc_cv.append(percent_cv)
    #add the percent cv column to the housekeeping datafram and output to file
    house_genes['Percent CV']=perc_cv
    house_genes.to_csv(sample_name+'_housekeeping.csv')
    delta_list=[]
    
    #looping through non-housekeeping samples
    for index,rows in array_samples.iterrows():
        x=rows['Sample Name']; y=rows['New Ct mean']; z=rows['Gene Symbol']
        appended_column=[] 
        #appended column will list the sample name, gene symbol, and later append the delta Ct
        appended_column.append(x)
        appended_column.append(z)
        
        #looping through housekeeping genes. Each sample will hit each housekeeping gene
        for index2, rows2 in house_genes.iterrows():
            f1=rows2['Sample Name']; f2=rows2['New Ct mean']
            #when a sample name from the non-housekeeping genes matches the housekeeping sample, it will calculate the delta Ct
            #and append to a list
            if x == f1:
                delta=y-f2
                appended_column.append(delta)
                
        # each loop will add a new row with sample name, sample name, and the delta Ct for each housekeeping gene
        delta_list.append(appended_column)
    #create a new pandas dataframe for the new list created when looping/matchin housekeeping genes and calculating delta Ct's
    delta_df=pandas.DataFrame(delta_list, columns=['Sample Name','Gene Symbol','B2M delta', 'GAPDH delta', 'HPRT1 delta', 'UBC assay1 delta', 
                                                   'UBC assay2 delta', 'YWHAZ delta'])
        
    #create a new column that calculates the mean for all the delat Ct's across housekeeping genes
    delta_df['Housekeeping mean']=delta_df[['B2M delta', 'GAPDH delta', 'HPRT1 delta', 'UBC assay1 delta', 
                                                   'UBC assay2 delta', 'YWHAZ delta']].mean(axis=1)
    
    #merge the delta Ct dataframe with the original sample input file
    output=pandas.merge(array_samples, delta_df, on=['Sample Name','Gene Symbol'], how='outer')
    # pulls all control sample names into new arrays
    inter_hek=output[output['Sample Name']=='HEK'] ; inter_uhr=output[output['Sample Name']=='UHR']
    # concats the two control sample arrays
    inter_df=pandas.concat([inter_hek,inter_uhr])
    # output with just the control samples for each sample
    inter_df.to_csv(sample_name+'_inter_controls.csv')
    
    output.to_csv(sample_name+'_output.csv')


# Takes in processed output files 
def array_qc(f1,f2,f3):
    ## read in files, change based on how many samples
    file1=pandas.read_csv(f1, usecols=[1,2,19]) 
    file2=pandas.read_csv(f2, usecols=[1,2,19]) 
    file3=pandas.read_csv(f3, usecols=[1,2,19]) 
    #file4=pandas.read_csv(f4, usecols=[1,2,19])
    #print(file1,file2,file3,file4)
    
    ## Change based on number of input files
    file_inputs=[file1,file2,file3]
    
    ## Change based on input file name
    samp1=f1.replace('_QuantStudio_export_processed_output.csv','') 
    samp2=f2.replace('_QuantStudio_export_processed_output.csv','') 
    samp3=f3.replace('_QuantStudio_export_processed_output.csv','')
    #samp4=f4.replace('_QuantStudio_export_processed_output.csv','')
    
    ## Change based on number of input files
    data=pandas.concat(file_inputs, keys=[samp1,samp2,samp3])
    
    #get all the gene names from the sample files
    genes=file1['Gene Symbol'].unique()
    HEK=[]
    UHR=[]
    # reset the index to the keys used when concat all input files
    data.reset_index(level=1, drop=True, inplace=True)
    
    # loop through the data and pull all HEK and UHR into separate arrays
    for j, k in data.iterrows():
        if k['Sample Name']=='HEK':
            HEK.append(k)
        elif k['Sample Name']=='UHR':
            UHR.append(k)
        else:
            pass
    # convert lists into dataframes
    hek_df=pandas.DataFrame(HEK); uhr_df=pandas.DataFrame(UHR)

    fig=plt.figure(figsize=(11,8))
    ax1= fig.add_subplot(111)
    
    # loop throught each grouped sample and get the list of all housekeeping means and plot on graph 1
    for i, x in hek_df.groupby(hek_df.index):
        graph_data=[]
        for k, v in x.iterrows():
            graph_data.append(v['Housekeeping mean'])
        ax1.plot(genes, graph_data, label=i, marker='o')
    
    plt.legend()
    plt.title('HEK control')
    plt.xticks(rotation=90); plt.xlabel('Genes'); plt.ylabel('Housekeeping mean')
    #plt.savefig('HEK control.png', dpi=300)
    
    # Repeated process for uhr identical to hek samples
    fig=plt.figure(figsize=(11,8))
    ax1= fig.add_subplot(111)
    for i, x in uhr_df.groupby(uhr_df.index):
        graph_data=[]
        for k, v in x.iterrows():
            graph_data.append(v['Housekeeping mean'])
        ax1.plot(genes, graph_data, label=i, marker='o')
        
    plt.legend()
    plt.title('UHR control')
    plt.xticks(rotation=90); plt.xlabel('Genes'); plt.ylabel('Housekeeping mean')
    plt.savefig('UHR control.png', dpi=300)
