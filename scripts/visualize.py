import pandas as pd
import matplotlib as mpl
# Use the PGF backend for LaTeX integration
mpl.use("pgf")
import matplotlib.pyplot as plt
import os
import numpy as np

# --- Configuration ---
# Define the pairs of benchmark files you want to compare.
# Each tuple contains: (path_to_diagonal_file, path_to_rowwise_file, name_for_output_plots)
BENCHMARK_PAIRS = [
    (
        "../benchmarks/results/diagonal_plain_matvec.txt",
        "../benchmarks/results/rowwise_plain_matvec.txt",
        "Plaintext Matrix-Vector Multiplication"
    ),
    (
        "../benchmarks/results/diagonal_enc_matvec.txt",
        "../benchmarks/results/rowwise_enc_matvec.txt",
        "Encrypted Matrix-Vector Multiplication"
    ),
    (
        "../benchmarks/results/diagonal_encmat_encmat.txt",
        "../benchmarks/results/rowwise_encmat_encmat.txt",
        "Encrypted Matrix-Matrix Multiplication"
    )
]

# Directory to save the generated plots and tables
OUTPUT_DIR = "../benchmarks/results/latex_visualizations"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Aesthetic Settings for LaTeX/PGF ---
# This setup configures matplotlib to use the system's LaTeX compiler
# and fonts for a seamless look in the final paper.
plt.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "font.family": "serif",
    "text.usetex": True,
    "pgf.rcfonts": False,
    "pgf.preamble": r"""
        \usepackage[utf8x]{inputenc}
        \usepackage[T1]{fontenc}
        \usepackage{booktabs}
    """,
    'font.size': 10,
    'axes.labelsize': 10,
    'legend.fontsize': 8,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
})

DIAG_COLOR = '#08519c'  # Deeper, professional blue
ROW_COLOR = '#d95f02'   # Professional orange

# --- Functions ---

def parse_benchmark_file(file_path):
    """
    Parses a benchmark text file into a pandas DataFrame.
    """
    if not os.path.exists(file_path):
        print(f"Warning: Benchmark file not found, skipping: {file_path}")
        return pd.DataFrame()

    data = []
    try:
        with open(file_path, "r") as f:
            lines = f.readlines()
            if len(lines) <= 1:
                return pd.DataFrame()

            for line in lines[1:]:
                parts = line.strip().split()
                if not parts:
                    continue
                
                first_numeric_idx = -1
                for i, part in enumerate(parts):
                    try:
                        float(part)
                        first_numeric_idx = i
                        break
                    except ValueError:
                        continue
                
                if first_numeric_idx != -1:
                    operation = " ".join(parts[:first_numeric_idx])
                    total_time = float(parts[first_numeric_idx])
                    calls = int(parts[first_numeric_idx + 1])
                    
                    data.append({
                        "Operation": operation,
                        "TotalTime(ms)": total_time,
                        "Calls": calls
                    })
    except Exception as e:
        print(f"Error parsing file {file_path}: {e}")
        return pd.DataFrame()
        
    return pd.DataFrame(data)

def generate_latex_barchart(df_compare, comparison_name):
    """
    Generates and saves a grouped bar chart as a .pgf file for LaTeX.
    """
    fig, ax = plt.subplots(figsize=(7.5, 4)) # Standard figure size for papers
    
    bar_width = 0.35
    x = np.arange(len(df_compare["Operation"]))

    ax.bar(x - bar_width/2, df_compare["TotalTime(ms)_diag"], width=bar_width, label="Diagonalwise", color=DIAG_COLOR)
    ax.bar(x + bar_width/2, df_compare["TotalTime(ms)_row"], width=bar_width, label="Rowwise", color=ROW_COLOR)

    ax.set_xticks(x)
    # Escape underscores for LaTeX
    labels = [op.replace('_', r'\_') for op in df_compare["Operation"]]
    ax.set_xticklabels(labels, rotation=45, ha="right")
    
    ax.set_ylabel("Total Time (ms)")
    ax.set_title(f"Total Time Comparison: {comparison_name}", weight='bold')
    ax.legend(frameon=False)
    ax.grid(axis='y', linestyle='--', color='#cccccc', alpha=0.7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    filename = os.path.join(OUTPUT_DIR, f"{comparison_name.replace(' ', '_')}_time_chart.pgf")
    plt.savefig(filename)
    plt.close()
    print(f"Generated LaTeX PGF chart: {filename}")

def generate_latex_table(df_compare, comparison_name):
    """
    Generates and saves a comparison table as a .tex file using booktabs style.
    """
    calls_table_df = df_compare[["Operation", "Calls_diag", "Calls_row"]]
    calls_table_df.columns = ["Operation", "Diagonalwise Calls", "Rowwise Calls"]
    calls_table_df = calls_table_df.astype({'Diagonalwise Calls': int, 'Rowwise Calls': int})

    # Escape special LaTeX characters in operation names
    calls_table_df["Operation"] = calls_table_df["Operation"].str.replace('_', r'\_', regex=False)

    # Begin LaTeX table string
    latex_string = r"""
\begin{table}[ht]
\centering
\caption{Operation Call Counts: """ + comparison_name + r"""}
\label{tab:""" + comparison_name.replace(' ', '_') + r"""}
\begin{tabular}{@{}lrr@{}}
\toprule
\textbf{Operation} & \textbf{Diagonalwise Calls} & \textbf{Rowwise Calls} \\
\midrule
"""
    # Add table rows
    for index, row in calls_table_df.iterrows():
        op, diag_calls, row_calls = row
        
        # Highlight the larger value with \bfseries
        diag_str = f"\\bfseries {diag_calls}" if diag_calls > row_calls else str(diag_calls)
        row_str = f"\\bfseries {row_calls}" if row_calls > diag_calls else str(row_calls)
        
        latex_string += f"{op} & {diag_str} & {row_str} \\\\\n"

    # End LaTeX table string
    latex_string += r"""\bottomrule
\end{tabular}
\end{table}
"""
    
    filename = os.path.join(OUTPUT_DIR, f"{comparison_name.replace(' ', '_')}_calls_table.tex")
    with open(filename, 'w') as f:
        f.write(latex_string)
    
    print(f"Generated LaTeX table: {filename}")


# --- Main Execution ---

if __name__ == "__main__":
    print(f"Starting benchmark visualization...")
    print(f"Output directory: {OUTPUT_DIR}")
    
    for diag_path, row_path, name in BENCHMARK_PAIRS:
        print("-" * 50)
        print(f"Processing pair: {name}")
        df_diag = parse_benchmark_file(diag_path)
        df_row = parse_benchmark_file(row_path)
        
        if not df_diag.empty and not df_row.empty:
            df_compare = pd.merge(
                df_diag,
                df_row,
                on="Operation",
                how="outer",
                suffixes=('_diag', '_row')
            ).fillna(0)
            
            generate_latex_barchart(df_compare, name)
            generate_latex_table(df_compare, name)
        else:
            print(f"Skipping visualization for '{name}' due to missing data.")

    print("-" * 50)
    print("Visualization script finished.")
    print("\nTo use the generated files in your LaTeX document:")
    print(r"1. Make sure you have `\usepackage{pgf}` and `\usepackage{booktabs}` in your preamble.")
    print(r"2. For charts, use: `\begin{figure}[ht] \centering \input{" + os.path.join(OUTPUT_DIR, "figure_name.pgf") + r"} \caption{Your caption} \end{figure}`")
    print(r"3. For tables, use: `\input{" + os.path.join(OUTPUT_DIR, "table_name.tex") + r"}`")

