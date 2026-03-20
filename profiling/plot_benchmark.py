import pandas as pd
import matplotlib.pyplot as plt
import os

# Configuration
# Configuration
DATA_FILE = 'profiling/benchmark_results_simdVSnosimd.csv'
OUTPUT_PLOT = 'report/figures/timing_plot.png'
SPEEDUP_PLOT = 'report/figures/speedup_plot.png'

def plot_benchmark():
    if not os.path.exists(DATA_FILE):
        print(f"Error: Data file {DATA_FILE} not found.")
        return

    try:
        # Read the CSV file
        df = pd.read_csv(DATA_FILE)
        
        # Separate data
        df_simd = df[df['simd'] == 1].sort_values('domain_size')
        df_nosimd = df[df['simd'] == 0].sort_values('domain_size')
        
        # Create the timing plot
        plt.figure(figsize=(10, 6))
        
        if not df_nosimd.empty:
            plt.plot(df_nosimd['domain_size'], df_nosimd['computation_time'], 'o-', linewidth=2, markersize=8, label='No SIMD')
        
        if not df_simd.empty:
            plt.plot(df_simd['domain_size'], df_simd['computation_time'], 's-', linewidth=2, markersize=8, label='SIMD')

        # Add labels and title
        plt.title('LBM Computation Time: SIMD vs No SIMD', fontsize=14)
        plt.xlabel('Domain Size (NxN)', fontsize=12)
        plt.ylabel('Computation Time (seconds)', fontsize=12)
        plt.legend()
        plt.grid(True, linestyle='--', alpha=0.7)
        
        # Save the timing plot
        plt.savefig(OUTPUT_PLOT)
        print(f"Timing plot saved to {OUTPUT_PLOT}")

        # Create Speedup Plot if we have both data
        if not df_simd.empty and not df_nosimd.empty:
             # Merge on domain_size to ensure we compare same sizes
             merged = pd.merge(df_nosimd, df_simd, on='domain_size', suffixes=('_nosimd', '_simd'))
             merged['speedup'] = merged['computation_time_nosimd'] / merged['computation_time_simd']
             
             plt.figure(figsize=(10, 6))
             plt.plot(merged['domain_size'], merged['speedup'], '^-', color='green', linewidth=2, markersize=8)
             plt.title('LBM SIMD Speedup Factor', fontsize=14)
             plt.xlabel('Domain Size (NxN)', fontsize=12)
             plt.ylabel('Speedup (Factor)', fontsize=12)
             plt.grid(True, linestyle='--', alpha=0.7)
             
             plt.savefig(SPEEDUP_PLOT)
             print(f"Speedup plot saved to {SPEEDUP_PLOT}")

        
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    plot_benchmark()
