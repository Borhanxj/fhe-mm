import pandas as pd
import matplotlib.pyplot as plt
import os

# ------------------- File paths -------------------
diag_path = "../benchmarks/results/diagonal_encmat_encmat.txt"
row_path = "../benchmarks/results/rowwise_encmat_encmat.txt"
output_dir = "../benchmarks/results"
os.makedirs(output_dir, exist_ok=True)

# ------------------- Optional operation rename mapping -------------------
rename_map = {
    # Row-wise
    "Encrypted matvec (column of B)": "MatVec (Row)",
    "Masking": "Mask (Row)",
    "Addition (row result)": "Add (Row)",
    "Rotate (row result)": "Rotate (Row)",
    "Addition (results)": "Add Final (Row)",
    "Mult + Relin + rescale": "Mult+Relin+Rescale (Row)",
    "Encode mask": "Encode Mask (Row)",
    "ModSwitch + Multiply Plain + Rescale": "ModSwitch+Mult+Rescale (Row)",
    
    # Diagonal-wise
    "diagonal matvec": "MatVec (Diag)",
    "Addition (result)": "Add (Diag)",
    "Rotation (enc vector)": "Rotate (Diag)",
    "Mult + Relin + Rescale": "Mult+Relin+Rescale (Diag)"
}

# ------------------- Parse benchmark file -------------------
def parse_benchmark(file_path):
    data = []
    with open(file_path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if not parts or parts[0] == "Operation":
                continue
            # Find where the numbers start
            i = 0
            while not parts[i].replace('.', '', 1).isdigit():
                i += 1
            operation = " ".join(parts[:i])
            total_time = float(parts[i])
            calls = int(parts[i + 1])
            avg_time = float(parts[i + 2])
            data.append({
                "Operation": rename_map.get(operation, operation),
                "TotalTime(ms)": total_time,
                "Calls": calls,
                "AvgTime(us)": avg_time
            })
    return pd.DataFrame(data)

# ------------------- Load data -------------------
df_diag = parse_benchmark(diag_path)
df_row = parse_benchmark(row_path)

# Normalize operation names for merging
def base_operation(op_name):
    return op_name.split()[0]  # first word: MatVec, Add, Rotate, etc.

df_diag["BaseOperation"] = df_diag["Operation"].apply(base_operation)
df_row["BaseOperation"] = df_row["Operation"].apply(base_operation)

# ------------------- Merge -------------------
df_compare = pd.merge(
    df_diag[["BaseOperation", "Calls", "TotalTime(ms)"]].rename(
        columns={"Calls": "Diagonalwise Calls", "TotalTime(ms)": "Diagonalwise Time"}
    ),
    df_row[["BaseOperation", "Calls", "TotalTime(ms)"]].rename(
        columns={"Calls": "Rowwise Calls", "TotalTime(ms)": "Rowwise Time"}
    ),
    on="BaseOperation",
    how="outer"
).fillna(0)

# ------------------- 1. Grouped Bar Chart for Total Time -------------------
plt.figure(figsize=(10, 6))
bar_width = 0.35
x = range(len(df_compare))

plt.bar([i - bar_width/2 for i in x], df_compare["Diagonalwise Time"],
        width=bar_width, label="Diagonalwise", color="#4C72B0")
plt.bar([i + bar_width/2 for i in x], df_compare["Rowwise Time"],
        width=bar_width, label="Rowwise", color="#55A868")

plt.xticks(ticks=x, labels=df_compare["BaseOperation"], rotation=30, ha="right")
plt.ylabel("Total Time (ms)")
plt.title("Total Time Comparison")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "total_time_comparison.png"))
plt.close()

# ------------------- 2. Calls Comparison Table -------------------
fig_height = max(2, 0.6 * len(df_compare))
fig, ax = plt.subplots(figsize=(10, fig_height))
ax.axis('off')

calls_table = df_compare[["BaseOperation", "Diagonalwise Calls", "Rowwise Calls"]]

tbl = ax.table(
    cellText=calls_table.values,
    colLabels=["Operation", "Diagonalwise Calls", "Rowwise Calls"],
    cellLoc='center',
    loc='center'
)

tbl.auto_set_font_size(False)
tbl.set_fontsize(10)
tbl.scale(1, 1.5)

# Highlight the higher value in red
for i in range(len(calls_table)):
    diag_calls = calls_table.iloc[i]["Diagonalwise Calls"]
    row_calls = calls_table.iloc[i]["Rowwise Calls"]
    if diag_calls > row_calls:
        tbl[(i + 1, 1)].set_facecolor("#ffcdd2")
        tbl[(i + 1, 2)].set_facecolor("white")
    elif row_calls > diag_calls:
        tbl[(i + 1, 2)].set_facecolor("#ffcdd2")
        tbl[(i + 1, 1)].set_facecolor("white")
    else:
        tbl[(i + 1, 1)].set_facecolor("white")
        tbl[(i + 1, 2)].set_facecolor("white")

plt.title("Operation Calls Comparison")
plt.subplots_adjust(left=0.2, right=0.8, top=0.9, bottom=0.1)
plt.savefig(os.path.join(output_dir, "calls_comparison.png"))
plt.close()
