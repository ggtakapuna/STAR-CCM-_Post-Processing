
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
        csv_1 = filedialog.askopenfilename(
            title="Select FIRST CSV file",
            filetypes=[("CSV files", "*.csv")],
            parent=root
        )
        if not csv_1:
            return

        csv_2 = filedialog.askopenfilename(
            title="Select SECOND CSV file",
            filetypes=[("CSV files", "*.csv")],
            parent=root
        )
        if not csv_2:
            return

        x1, cp1 = load_csv(csv_1)
        x2, cp2 = load_csv(csv_2)

    except (OSError, ValueError) as error:
        messagebox.showerror("CSV Error", str(error), parent=root)
        return
    finally:
        root.destroy()

    fig, ax = plt.subplots(figsize=(11, 6))

    ax.scatter(x1, cp1, s=12, alpha=0.7,
               label=Path(csv_1).stem)
    ax.scatter(x2, cp2, s=12, alpha=0.7,
               label=Path(csv_2).stem)

    ax.set_title("Pressure Coefficient Comparison", fontsize=16)
    ax.set_xlabel("Position X (m)", fontsize=12)
    ax.set_ylabel("Pressure Coefficient, Cp", fontsize=12)

    ax.grid(True, linestyle="--", alpha=0.35)
    ax.legend(fontsize=8)

    # Uncomment for conventional airfoil Cp plots:
    # ax.invert_yaxis()

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
