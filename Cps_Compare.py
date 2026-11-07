
import csv
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt


def load_csv(path):
    x_values = []
    cp_values = []

    with open(path, "r", newline="", encoding="utf-8-sig") as file:
        reader = csv.reader(file)
        next(reader, None)  # Skip column headings

        for row in reader:
            if not row or all(not item.strip() for item in row):
                continue

            x_values.append(float(row[0]))
            cp_values.append(float(row[1]))

    return x_values, cp_values


def main():
    root = tk.Tk()
    root.withdraw()

    try:
        # Select multiple CSV files
        csv_files = filedialog.askopenfilenames(
            title="Select Cp CSV Files",
            filetypes=[("CSV files", "*.csv")],
            parent=root
        )

        if not csv_files:
            return

        plot_count = 0

        for csv_file in csv_files:

            try:
                x, cp = load_csv(csv_file)

                if not x:
                    print(f"No data found in {csv_file}")
                    continue

                # Create a separate figure for each CSV
                fig, ax = plt.subplots(
                    figsize=(11, 6),
                    num=Path(csv_file).stem
                )

                ax.scatter(
                    x, cp,
                    s=12,
                    alpha=0.7,
                    color="tab:blue",
                    label=Path(csv_file).stem
                )

                ax.set_title(
                    Path(csv_file).stem,
                    fontsize=16
                )

                ax.set_xlabel(
                    "Position X (m)",
                    fontsize=12
                )

                ax.set_ylabel(
                    "Pressure Coefficient, Cp",
                    fontsize=12
                )

                ax.grid(
                    True,
                    linestyle="--",
                    alpha=0.35
                )

                ax.legend(fontsize=8)

                # Uncomment for conventional airfoil Cp plots:
                # ax.invert_yaxis()

                fig.tight_layout()

                plot_count += 1

                print(f"Loaded: {Path(csv_file).name}")

            except (OSError, ValueError, IndexError) as error:
                print(f"Error loading {csv_file}: {error}")

        print(f"\nSuccessfully created {plot_count} Cp plots.")

        if plot_count > 0:
            # Display all figures
            plt.show()
        else:
            messagebox.showerror(
                "No Plots",
                "No valid CSV files were found.",
                parent=root
            )

    finally:
        root.destroy()


if __name__ == "__main__":
    main()
