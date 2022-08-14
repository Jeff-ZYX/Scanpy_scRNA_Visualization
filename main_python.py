# %% [markdown]
# Written by Jeff Zhang

# %%
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.tools as tl 
import os

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

# %%

def show_dup(data):
    '''
    Prints out number of duplicate genes and gene names in data file
    '''
    print("\n----- Showing Duplicates -----\n")
    dup_genes = data.var_names[data.var_names.duplicated()].tolist()
    print(dup_genes)
    print(len(dup_genes))

# %%
def load_file(path,type=""):
    '''
    Loads data files into Anndata format
    For mtx data 
      - 3 files are required with the EXACT following filenames("matrix.mtx","barcodes.tsv","genes.tsv") for the function to run. 
      - All 3 files must be in one folder which the "path" variable refers to.
    For h5 data
      - The "path" variable refers to the h5 file.
    '''
    if type == "":
        filetype = input("What type of file are you loading? \nThis program can only read h5 and mtx files \nIf reading mtx files, Ensure the file_path provided is a folder that contains a 'barcodes.tsv', 'genes.tsv', and 'matrix.mtx' file")
    print("\n----- Loading File -----\n")
    if type == "mtx":
        adata = sc.read_10x_mtx(path,var_names='gene_symbols')
        show_dup(adata)
        
        return(adata)
    elif type == "h5":
        adata = sc.read_10x_h5(path)
        show_dup(adata)
        
        return(adata)
    
    


# %%
def min_filter(adata,gene_n=False,cell_n=False):
    '''
    Filters the anndata input and returns filtered data
       - exclude cells with less than "gene_n" number of genes
       - exclude genes which appear in less than "cell_n" number of cells
       - calculate the percentage of mitocondrial genes per cell and add to the metadata
    '''
    print("\n----- Filtering by minimum gene and cell count required -----\n")
    if gene_n != False:
        
        sc.pp.filter_cells(adata, min_genes=gene_n)
    elif cell_n != False:
        
        sc.pp.filter_genes(adata, min_cells=cell_n)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    return(adata)

    
    
    


# %%
def max_filter(adata,gene_n=False,mito_n=False):
    '''
    Filters out cells above certain cut offs and returns filtered data
        - filters out cells with more thant "gene_n" number of genes
        - filters out cells with less than "mito_n"% mito gene umi count
    '''
    print("\n----- Filtering by gene and mitochondria count cutoff -----\n")
    if gene_n != False:
        adata = adata[adata.obs.n_genes_by_counts < gene_n, :] 
    if mito_n != False:
        adata = adata[adata.obs.pct_counts_mt < mito_n, :]
    return adata

# %%
def hvg(adata, minM, maxM, minD):
    '''
    DEFUNCT
    Differentiates and Isolates Highly Variable Genes from adata and returns filtered adata
    
    '''
    print("\n----- Differentiating Highly Variable Genes -----\n")
    sc.pp.highly_variable_genes(adata, min_mean=minM, max_mean=maxM, min_disp = minD)
    print("\n----- Highly Variable Genes Count -----\n")
    print(adata.var.highly_variable.value_counts())
    print("\n----- Printing chart of Highly/Non-Highly Variable Genes -----\n")
    sc.pl.highly_variable_genes(adata)
    print("\n----- Isolating Highly Variable Genes -----\n")
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    return adata


# %%
def plot_principal(adata,prin_comp):
    '''
    DEFUNCT
    For every principal component in "prin_comp", this function plots it
    '''
    print("\n----- Plotting Principal components -----\n")
    for i in prin_comp:
        try:
            sc.pl.pca(adata,color=str(i))
        except:
            print(f"an error occurred with {i}")
            pass

# %%
def main(path,full_list,type,hvg_min_mean=0.0125,hvg_max_mean=3,hvg_min_D=0.5,top=20,gene_num=False,cell_num=False,max_gene=False,max_mito=False,variance_max=10,prin_gene=[],n_n=None,n_p=None):
    '''
    A reproducible function for the preprocessing of scRNA files followed by visualization for 
    genes of interest

    -Explanation of variables tbc...
    
    '''
    print(f"one: {full_list}")
    adata = load_file(path,type)
    adata.var_names_make_unique()
    
    #Plot "top" number of expressed genes
    sc.pl.highest_expr_genes(adata,n_top=top)

    adata = min_filter(adata,gene_num,cell_num)

    #Plot violin and scatter plots for genes
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

    adata = max_filter(adata=adata,gene_n=max_gene,mito_n=max_mito)

    #normalize and scale data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    #filter out highly variable genes and plot them
    sc.pp.highly_variable_genes(adata, min_mean=hvg_min_mean, max_mean=hvg_max_mean, min_disp = hvg_min_D)
    sc.pl.highly_variable_genes(adata)

    #Include and filter out necessary data for plotting
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]

    #Regress and scale data
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)

    #Compute PCA Matrix
    sc.tl.pca(adata, svd_solver='arpack')

    #Plot PCA for principal genes
    for i in prin_gene:
        try:
            sc.pl.pca(adata,color=str(i))
        except:
            print(f"an error occurred with {i}")
            pass
    
    #Plot variance ratio
    sc.pl.pca_variance_ratio(adata, log=True)

    #Find Neighbours for data points
    sc.pp.neighbors(adata, n_neighbors=n_n, n_pcs=n_p)

    #Employ Leiden Algorithm on Data
    sc.tl.leiden(adata)

    # Mapping out the coarse-grained connectivity structures of complex manifolds using partition-based graph abstraction (PAGA)
    tl.paga(adata) 
    sc.pl.paga(adata, plot=False)
    tl.umap(adata, init_pos='paga')

    #Plot UMAP of principal Gene
    for i in prin_gene:
        try:
            sc.pl.umap(adata,color=str(i))
        except:
            print(f"an error occurred with {i}")
            pass
    
    full_list_2 = full_list.copy()
    full_list_2.append("leiden")

    #Plot UMAP and violin chart of Genes in question to visualize relative expression
    sc.pl.umap(adata, color=full_list_2, wspace=0.4, color_map="Blues")
    sc.pl.violin(adata, full_list, groupby='leiden')
    return adata
    

    

# %% [markdown]
# 3 Examples of data processing visualization from datasets in the following github repository:
# https://github.com/Jeff-ZYX/Scanpy_scRNA_Visualization.git
# 
# The Github repository aims to replicate R processing of data in the following paper with Scanpy:
# https://link.springer.com/content/pdf/10.1007/s11684-020-0754-0.pdf
# 

# %%
#Lung data
file_type ="h5"
file_path = ".\data\Lung\GSE122960_RAW\GSM3489182_Donor_01_filtered_gene_bc_matrices_h5.h5"
gene_min = 200
cell_min = 3
gene_cutoff = 3000
mito_cutoff = 10
hvg_minMean=0.01,
hvg_maxMean=8,
hvg_minD=0.5
principal_components = ['SFTPB']
num_n = 20
num_p = 30
interest = ['SFTPB','SFTPC','ACE2']


adata = main(path=file_path,type=file_type,gene_num=gene_min,cell_num=cell_min,max_gene=gene_cutoff,max_mito=mito_cutoff,prin_gene=principal_components,n_n=num_n,n_p=num_p,full_list=interest,hvg_min_mean=hvg_minMean,hvg_max_mean=hvg_maxMean,hvg_min_D=hvg_minD)


# %%
# #Bladder Data
file_type = "mtx"
file_path = ".\data\Bladder\collated_data"
gene_min = 200
cell_min = 3
gene_cutoff = 2500
mito_cutoff = 5
principal_components = ['ACE2']
num_n = 10
num_p = 40
interest = ['SPINK1','CLDN4','ACE2']

main(path=file_path,type=file_type,gene_num=gene_min,cell_num=cell_min,max_gene=gene_cutoff,max_mito=mito_cutoff,prin_gene=principal_components,n_n=num_n,n_p=num_p,full_list=interest)

# %%
#Ileum data
file_type = "mtx"
file_path = ".\data\Ileum\collated_data"
gene_min = 200
cell_min = 3
gene_cutoff = 2500
mito_cutoff = 5
principal_components = ['CST3']
num_n = 20
num_p = 30
interest = ['FABP6','ACE2','ANPEP']

main(path=file_path,type=file_type,gene_num=gene_min,cell_num=cell_min,max_gene=gene_cutoff,max_mito=mito_cutoff,prin_gene=principal_components,n_n=num_n,n_p=num_p,full_list=interest)


