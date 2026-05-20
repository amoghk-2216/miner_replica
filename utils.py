import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from lifelines import CoxPHFitter
from scipy.stats import binomtest
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
import matplotlib.pyplot as plt


def validate_coexpression(cluster, expr, threshold=0.30, min_genes = 2):                    
    
    '''

    Args:           cluster (list)      : list of gene names
                    expr (df)           : df containing seq data (genes * patients)
                    threshold (float)   : % of variance that must be explained by PC1; default is 30%
                    min_genes (int)     : minimum number of genes required for PCA (default 2)


    Returns:        bool                : True if variance explained >= threshold


    '''
    
    try:
        
        valid_genes = [gene for gene in cluster if gene in expr.index]       # Only considers genes in data fed to function
        num_genes = len(valid_genes)
    
        # PCA requires at least 2 features (genes) to evaluate co-expression
        if num_genes < 2:
            return False

        try:
            # Extract submatrix. Since expr is (genes, patients), 
            # submatrix.T results in shape: (patients, genes)
            submatrix = expr.loc[valid_genes].T 
        
            # Ensure we have at least 2 patients
            if submatrix.shape[0] < 2:
                return False
                
            # Fit 1-component PCA
            pca = PCA(n_components=1, random_state=42)
            pca.fit(submatrix)
            variance_explained = pca.explained_variance_ratio_[0]           # returns a list and extracts PC1 which is its first value

            return variance_explained >= threshold                          # does PC1 cause >= threshold variance?
        
        except Exception:
        # Handles mathematical edge cases
            return False
            
        
    except Exception as e:
        return False




def get_eigengene(regulon, expr):

    '''
    
    Args:           regulon (list)      : list of gene names in regulon          
                    expr (df)           : dataframe containing seq data (genes * patients)

    Returns:        eigen (series)      : eigengene score for a particular regulon indexed by patient id

    '''
    
    try:
        
        # Filter down to the genes that are present in the expression matrix 
        # This is redundant here since validated only contains genes in expr
        valid_genes = [gene for gene in regulon if gene in expr.index]
        
        if len(valid_genes) < 2:                          # PCA needs >= 2 genes
            return None

        submatrix = expr.loc[valid_genes].T               # Extracts cluster genes, transposes shape to patients x genes as this is the shape PCA needs
                                                          # This is because we need one score per patient i.e patients are rows/ records/ observations

        pca = PCA(n_components=1, random_state=42)               
        scores = pca.fit_transform(submatrix)             # Fits PCA and get PC1 score for each patient
        eigen = pd.Series(                                # Convert to flattened Series with patient IDs as index
            scores.flatten(),
            index=submatrix.index                  
        )
        
        mean_expr = submatrix.mean(axis=1)                # Mean expression of cluster genes per patient
        if np.corrcoef(eigen, mean_expr)[0,1] < 0:        # .corrcoeff generates n*n matrix where offdiagonal elements are correlation coeffs
                                                          # Here correlation is symmetric i.e matrix[0, 1] = matrix[1,0] = corrcoef  
            
            eigen = -eigen                                # Flip sign (so that high exp relative to mean gives eigengene > 0)
            
        return eigen                                      # Returns eigengene value for a particular regulon
    
    except Exception as e:
        return None                                       # Return None if anything goes wrong




def cox_regulon(eigen, clinical, duration_col='OS_MONTHS', event_col='OS_STATUS'):    

    '''
    
    Args :          eigen (pd. Series)  :   output of get_eigengene
                    clinical (df)       :   dataframe containing clinical data
                    duration_col (str)  :   column name for survival duration
                    event_col (str)     :   column name for survival status

    Returns :       hr (float)          :   hazard risk ratio ( > 1 indicates high-risk) 
                    pvalue (float)      :   p-value for Cox regression
    
    '''
     
    try:
        df = pd.DataFrame({                                                           # 3-column df w/ eigengene, survival duration and survival status
            'eigen': eigen, 
            'duration': clinical[duration_col],
            'event': clinical[event_col]
        }).dropna()
        
        if len(df) < 10:                                                              # if < 10 non-NaN patients, not enough data
            return None, None
        
        cph = CoxPHFitter()                                                           # model basically asks if eigengene predicts how quickly patient dies
        cph.fit(df, duration_col='duration', event_col='event')
        
        hr = cph.hazard_ratios_['eigen']
        pval = cph.summary['p']['eigen']
        
        return hr, pval                                                               # returns hazard ratio and p-value
    
    except Exception as e:
        return None, None
    



def quantize_matrix_fast(expr_df, disease_regulons, validated):

    '''

    Args:           expr_df (df)                    : expression matrix of shape (genes x patients), indexed by gene name, columns are patient IDs
                    disease_regulons (list[int])    : list of regulon indices to quantize                                
                    validated (list[list])          : list of lists where validated[i] is the list of genes in regulon i
                                          
    Returns:        quantized_matrix (df)           : (regulons x patients), indexed by regulon index, columns are patient IDs
                                                      values are:
                                                      1  → regulon significantly overexpressed in patient
                                                      -1  → regulon significantly underexpressed in patient
                                                      0  → no significant activity

    '''
    
    patients = expr_df.columns.tolist()
    n_genes = len(expr_df)
    
    # Precompute gene index lookup
    gene_to_idx = {g: i for i, g in enumerate(expr_df.index)}
    
    # Precompute regulon gene indices
    regulon_gene_idx = {}
    for reg_idx in disease_regulons:
        genes = validated[reg_idx]
        idx_list = [gene_to_idx[g] for g in genes if g in gene_to_idx]
        if len(idx_list) >= 3:
            regulon_gene_idx[reg_idx] = np.array(idx_list)
    
    valid_regulons = list(regulon_gene_idx.keys())
    n_regulons = len(valid_regulons)
    n_patients = len(patients)
    
    result = np.zeros((n_regulons, n_patients), dtype=np.int8)
    
    # Precompute thresholds
    lower_thresh = n_genes / 3
    upper_thresh = 2 * n_genes / 3
    
    expr_values = expr_df.values  # genes × patients numpy array
    
    for j in range(n_patients):
        if j % 50 == 0:
            print(f"Patient {j}/{n_patients}...")
        
        # Rank all genes for this patient in one shot
        ranks = pd.Series(expr_values[:, j]).rank(method='average').values
        
        for i, reg_idx in enumerate(valid_regulons):
            gene_idx = regulon_gene_idx[reg_idx]
            reg_ranks = ranks[gene_idx]
            n = len(reg_ranks)
            
            upper = int((reg_ranks > upper_thresh).sum())
            lower = int((reg_ranks < lower_thresh).sum())
            
            p_upper = binomtest(upper, n, 1/3, alternative='greater').pvalue
            p_lower = binomtest(lower, n, 1/3, alternative='greater').pvalue
            
            if p_upper <= 0.05:
                result[i, j] = 1
            elif p_lower <= 0.05:
                result[i, j] = -1
    
    quantized_matrix = pd.DataFrame(result, index=valid_regulons, columns=patients)
    return quantized_matrix




def compute_guan_rank(survival_df, duration_col='OS_MONTHS', event_col='OS_STATUS'):
    '''
    Args:
        survival_df  (pd.DataFrame) : indexed by PATIENT_ID with duration and event cols
        duration_col (str)          : survival time column
        event_col    (str)          : event indicator (1=dead, 0=censored)

    Returns:
        guan_scores  (pd.Series)    : GuanRank scores indexed by PATIENT_ID
    '''
    df = survival_df[[duration_col, event_col]].copy().reset_index()
    n = len(df)
    scores = np.zeros(n)

    for i in range(n):
        t_i = df[duration_col].iloc[i]
        e_i = df[event_col].iloc[i]
        if e_i == 1:
            scores[i] = (df[duration_col] > t_i).sum() / (n - 1)
        else:
            scores[i] = 1 - ((df[duration_col] < t_i) & (df[event_col] == 1)).sum() / (n - 1)

    df['GuanRank'] = scores
    return df.set_index(df.columns[0])['GuanRank']




def km_plot(risk_scores, survival_df, title='KM curve',
            duration_col='OS_MONTHS', event_col='OS_STATUS', save_path=None):
    '''
    Args:
        risk_scores  (pd.Series)    : risk scores indexed by PATIENT_ID
        survival_df  (pd.DataFrame) : survival data indexed by PATIENT_ID
        title        (str)          : plot title
        duration_col (str)          : survival time column
        event_col    (str)          : event indicator column
        save_path    (str)          : path to save figure (optional)

    Returns:
        pval         (float)        : log-rank p-value
    '''
    high_risk = risk_scores[risk_scores > risk_scores.median()].index
    low_risk  = risk_scores[risk_scores <= risk_scores.median()].index

    kmf = KaplanMeierFitter()
    fig, ax = plt.subplots(figsize=(8, 6))

    kmf.fit(survival_df.loc[high_risk, duration_col],
            survival_df.loc[high_risk, event_col],
            label=f'High risk (n={len(high_risk)})')
    kmf.plot_survival_function(ax=ax, color='red')

    kmf.fit(survival_df.loc[low_risk, duration_col],
            survival_df.loc[low_risk, event_col],
            label=f'Low risk (n={len(low_risk)})')
    kmf.plot_survival_function(ax=ax, color='blue')

    results = logrank_test(
        survival_df.loc[high_risk, duration_col],
        survival_df.loc[low_risk,  duration_col],
        event_observed_A=survival_df.loc[high_risk, event_col],
        event_observed_B=survival_df.loc[low_risk,  event_col]
    )

    ax.set_title(f'{title} (p={results.p_value:.4f})')
    ax.set_xlabel('Time (months)')
    ax.set_ylabel('Survival probability')
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150)

    plt.show()
    print(f"Log-rank p-value: {results.p_value:.4e}")
    return results.p_value




def concordance(risk_scores, survival_df,
                duration_col='OS_MONTHS', event_col='OS_STATUS'):
    '''
    Args:
        risk_scores  (pd.Series)    : risk scores indexed by PATIENT_ID
                                      higher score = higher risk = shorter survival
        survival_df  (pd.DataFrame) : survival data indexed by PATIENT_ID
        duration_col (str)          : survival time column
        event_col    (str)          : event indicator column

    Returns:
        ci           (float)        : concordance index
    '''
    common = risk_scores.index.intersection(survival_df.index)
    return concordance_index(
        survival_df.loc[common, duration_col],
        -risk_scores.loc[common],    # negate internally — higher risk = shorter survival
        survival_df.loc[common, event_col]
    )

