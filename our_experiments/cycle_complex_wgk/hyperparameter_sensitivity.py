import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np

# Get the script directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Set up the results directory (using absolute path)
RESULTS_DIR = os.path.join(SCRIPT_DIR, 'topowgk_outputs', 'processed_svm_results')

# Set output directory to the desired location
OUTPUT_DIR = os.path.join(SCRIPT_DIR, 'topowgk_outputs', 'visualizations')

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Define color and marker styles for each metric
METRICS = {
    'acc': {'color': 'blue', 'marker': 'o', 'label': 'ACC'},
    'f1': {'color': 'red', 'marker': 's', 'label': 'F1'},
    'roc_auc': {'color': 'green', 'marker': '^', 'label': 'ROC_AUC'}
}

# Define column mappings
COLUMN_MAPPING = {
    'acc': 'acc',
    'f1': 'f1_score',
    'roc_auc': 'roc_auc_score'
}

# Define datasets for appendix
APPENDIX_DATASETS = {
    'MSRC_9', 'ENZYMES', 'PROTEINS', 'COX2', 'BZR', 'DHFR_MD', 'NCI1', 'IMDB-BINARY'
}

def parse_filename(filename):
    """Parse filename to extract hyperparameter type and dataset name."""
    if filename.startswith('profile_alpha_'):
        param_type = 'alpha'
        dataset = filename.replace('profile_alpha_', '').replace('_svm_results.csv', '')
    elif filename.startswith('profile_gk_gamma_'):
        param_type = 'gk_gamma'
        dataset = filename.replace('profile_gk_gamma_', '').replace('_svm_results.csv', '')
    elif filename.startswith('profile_n_csg_layers_'):
        param_type = 'n_csg_layers'
        dataset = filename.replace('profile_n_csg_layers_', '').replace('_svm_results.csv', '')
    else:
        return None, None
    return param_type, dataset

def load_and_organize_data():
    """Load all profile files and organize by hyperparameter type and dataset."""
    data = {
        'alpha': {},
        'gk_gamma': {},
        'n_csg_layers': {}
    }
    
    # Get all profile files
    profile_files = [f for f in os.listdir(RESULTS_DIR) if f.startswith('profile_') and f.endswith('.csv')]
    
    for filename in profile_files:
        param_type, dataset = parse_filename(filename)
        if param_type is None:
            continue
            
        filepath = os.path.join(RESULTS_DIR, filename)
        df = pd.read_csv(filepath)
        
        # Determine which column contains the varying hyperparameter
        if param_type == 'alpha':
            param_col = 'alpha'
        elif param_type == 'gk_gamma':
            param_col = 'gk_gamma'
        elif param_type == 'n_csg_layers':
            param_col = 'n_csg_layers'
        
        # Group by the parameter column and calculate mean for each metric
        # This handles cases where multiple rows have the same parameter value
        agg_dict = {
            'acc': 'mean',
            'acc_std': 'mean',
            'f1_score': 'mean',
            'f1_score_std': 'mean',
            'roc_auc_score': 'mean',
            'roc_auc_score_std': 'mean'
        }
        
        df_aggregated = df.groupby(param_col, as_index=False).agg(agg_dict)
        
        # Sort by the parameter column to ensure proper line plotting
        df_aggregated = df_aggregated.sort_values(param_col).reset_index(drop=True)
        
        data[param_type][dataset] = df_aggregated
    
    return data

def plot_hyperparameter_sensitivity(data, param_type, datasets=None, is_appendix=False):
    """Create a grid plot for a specific hyperparameter type with specified datasets."""
    
    # Mapping from param_type to LaTeX symbols
    PARAM_SYMBOLS = {
        'alpha': r'$\beta$',
        'gk_gamma': r'$\gamma$',
        'n_csg_layers': r'$N_{hcc}$'
    }
    
    # Get list of datasets to plot
    if datasets is None:
        datasets = list(data[param_type].keys())
    
    # Filter datasets that exist in data
    datasets = [ds for ds in datasets if ds in data[param_type]]
    
    if not datasets:
        print(f"No datasets to plot for {param_type}")
        return
    
    # Calculate grid dimensions
    n_plots = len(datasets)
    n_cols = min(4, n_plots)
    n_rows = (n_plots + n_cols - 1) // n_cols
    
    # Create figure with subplots
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 5 * n_rows))
    if n_rows == 1 and n_cols == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    # Determine which column contains the varying hyperparameter
    if param_type == 'alpha':
        param_col = 'alpha'
    elif param_type == 'gk_gamma':
        param_col = 'gk_gamma'
    elif param_type == 'n_csg_layers':
        param_col = 'n_csg_layers'
    
    # Plot each dataset
    for i, dataset in enumerate(datasets):
        ax = axes[i]
        df = data[param_type][dataset]
        
        # Plot each metric
        for metric_key, metric_info in METRICS.items():
            col_name = COLUMN_MAPPING[metric_key]
            
            ax.plot(df[param_col], df[col_name], 
                    color=metric_info['color'], 
                    marker=metric_info['marker'], 
                    markersize=6, 
                    linestyle='-', 
                    linewidth=1.5,
                    label=metric_info['label'])
        
        ax.set_xlabel(PARAM_SYMBOLS[param_type], fontsize=15)
        ax.set_ylabel('Score', fontsize=15)
        ax.set_title(dataset, fontsize=15, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Use scientific notation for small numbers if needed
        if param_type == 'gk_gamma':
            ax.ticklabel_format(style='scientific', axis='x', scilimits=(0, 0))
    
    # Hide unused subplots
    for i in range(n_plots, len(axes)):
        axes[i].set_visible(False)
    
    # Add a unified legend for the entire figure
    handles = []
    labels = []
    for metric_info in METRICS.values():
        handles.append(mlines.Line2D([0], [0], color=metric_info['color'], 
                                     marker=metric_info['marker'], linestyle='-',
                                     markersize=10))
        labels.append(metric_info['label'])
    
    fig.legend(handles, labels, loc='upper center', 
              bbox_to_anchor=(0.5, 0.92), ncol=3, fontsize=15, frameon=True)
    
    # Set main title
    fig.suptitle(f'Hyperparameter Sensitivity: {PARAM_SYMBOLS[param_type]}', fontsize=20, fontweight='bold', y=0.95)
    
    plt.tight_layout(rect=(0, 0, 1, 0.92))
    
    # Determine file prefix
    prefix = 'supp_' if is_appendix else ''
    
    # Save the figure
    output_path = os.path.join(OUTPUT_DIR, f'{prefix}hyperparameter_sensitivity_{param_type}.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f'Saved: {output_path}')
    
    plt.close()

def save_data_to_csv(data, param_type, datasets=None, is_appendix=False):
    """Save aggregated data for a specific hyperparameter type to CSV."""
    
    # Determine which column contains the varying hyperparameter
    if param_type == 'alpha':
        param_col = 'alpha'
    elif param_type == 'gk_gamma':
        param_col = 'gk_gamma'
    elif param_type == 'n_csg_layers':
        param_col = 'n_csg_layers'
    
    # Get list of datasets to save
    if datasets is None:
        datasets = list(data[param_type].keys())
    
    # Filter datasets that exist in data
    datasets = [ds for ds in datasets if ds in data[param_type]]
    
    if not datasets:
        print(f"No datasets to save for {param_type}")
        return
    
    # Combine selected datasets into one DataFrame
    all_data = []
    
    for dataset in datasets:
        df = data[param_type][dataset]
        # Make a copy and add dataset column
        df_copy = df.copy()
        df_copy['dataset'] = dataset
        
        # Select and reorder columns
        cols = ['dataset', param_col, 'acc', 'f1_score', 'roc_auc_score']
        df_copy = df_copy[cols]
        
        all_data.append(df_copy)
    
    # Concatenate all datasets
    combined_df = pd.concat(all_data, ignore_index=True)
    
    # Determine file prefix
    prefix = 'supp_' if is_appendix else ''
    
    # Save to CSV
    output_path = os.path.join(OUTPUT_DIR, f'{prefix}hyperparameter_sensitivity_{param_type}_data.csv')
    combined_df.to_csv(output_path, index=False)
    print(f'Saved data: {output_path}')

def main():
    print("Loading data...")
    data = load_and_organize_data()
    
    # Print summary
    for param_type in ['alpha', 'gk_gamma', 'n_csg_layers']:
        print(f"\n{param_type.upper()}: {len(data[param_type])} datasets")
        print(f"  Datasets: {', '.join(sorted(data[param_type].keys()))}")
    
    print("\nSaving data to CSV...")
    
    # Save data to CSV for each hyperparameter type and section
    for param_type in ['alpha', 'gk_gamma', 'n_csg_layers']:
        all_datasets = list(data[param_type].keys())
        main_datasets = [ds for ds in all_datasets if ds not in APPENDIX_DATASETS]
        supp_datasets = [ds for ds in all_datasets if ds in APPENDIX_DATASETS]
        
        if main_datasets:
            print(f"\nSaving {param_type} data (main)...")
            save_data_to_csv(data, param_type, datasets=main_datasets, is_appendix=False)
        
        if supp_datasets:
            print(f"\nSaving {param_type} data (appendix)...")
            save_data_to_csv(data, param_type, datasets=supp_datasets, is_appendix=True)
    
    print("\nGenerating plots...")
    
    # Generate plots for each hyperparameter type and section
    for param_type in ['alpha', 'gk_gamma', 'n_csg_layers']:
        all_datasets = list(data[param_type].keys())
        main_datasets = [ds for ds in all_datasets if ds not in APPENDIX_DATASETS]
        supp_datasets = [ds for ds in all_datasets if ds in APPENDIX_DATASETS]
        
        if main_datasets:
            print(f"\nProcessing {param_type} (main)...")
            plot_hyperparameter_sensitivity(data, param_type, datasets=main_datasets, is_appendix=False)
        
        if supp_datasets:
            print(f"\nProcessing {param_type} (appendix)...")
            plot_hyperparameter_sensitivity(data, param_type, datasets=supp_datasets, is_appendix=True)
    
    print("\nAll plots and CSV files generated successfully!")
    print(f"Output directory: {OUTPUT_DIR}")

if __name__ == '__main__':
    main()
