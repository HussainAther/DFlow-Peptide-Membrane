import matplotlib.pyplot as plt
import pandas as pd
import json

class Visualization:
    @staticmethod
    def plot_membrane_thickness(filename="data/membrane_thickness.json"):
        """Plot the evolution of membrane thickness over time."""
        with open(filename, "r") as f:
            data = json.load(f)
        
        plt.figure(figsize=(8, 5))
        plt.plot(range(len(data)), [data["Membrane_Thickness"]] * len(data), label="Membrane Thickness", color='red')
        plt.xlabel("Cycles")
        plt.ylabel("Membrane Thickness")
        plt.title("Membrane Thickness Evolution Over Time")
        plt.legend()
        plt.show()
    
    @staticmethod
    def plot_peptide_selection(filename="data/peptide_evolution.csv"):
        """Plot the number of peptides surviving each cycle."""
        df = pd.read_csv(filename)
        
        plt.figure(figsize=(8, 5))
        plt.plot(df["Peptide_Count"], label="Surviving Peptides", color='blue')
        plt.xlabel("Cycles")
        plt.ylabel("Number of Peptides")
        plt.title("Peptide Selection Over Time")
        plt.legend()
        plt.show()

if __name__ == "__main__":
    Visualization.plot_membrane_thickness()
    Visualization.plot_peptide_selection()

