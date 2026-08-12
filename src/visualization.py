import matplotlib.pyplot as plt
import pandas as pd
import json

class Visualization:
    @staticmethod
    def plot_membrane_thickness(save_path="membrane_thickness.png"):
        df = pd.read_csv("data/peptide_evolution.csv")
        plt.figure(figsize=(8, 5))
        plt.plot(df["Cycle"], df["Membrane_Thickness"], label="Membrane Thickness")
        plt.xlabel("Cycle")
        plt.ylabel("Thickness")
        plt.title("Membrane Thickness Over Time")
        plt.legend()
        plt.tight_layout()
        plt.savefig(save_path, dpi=300)
        plt.close()   
 
    @staticmethod
    def plot_peptide_selection(save_path="peptide_selection.png"):
        df = pd.read_csv("data/peptide_evolution.csv")
        plt.figure(figsize=(8, 5))
        plt.plot(df["Cycle"], df["Peptide_Count"], label="Peptide Count", color='orange')
        plt.xlabel("Cycle")
        plt.ylabel("Count")
        plt.title("Peptide Selection Over Time")
        plt.legend()
        plt.tight_layout()
        plt.savefig(save_path, dpi=300)
        plt.close()

if __name__ == "__main__":
    Visualization.plot_membrane_thickness()
    Visualization.plot_peptide_selection()

