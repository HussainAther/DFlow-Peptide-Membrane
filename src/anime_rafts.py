import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from pathlib import Path

df = pd.read_csv("logs_v2/vesicle_state_log.csv")
outpath = Path("logs_v2/plots/raft_growth.gif")

fig, ax = plt.subplots()
line1, = ax.plot([], [], label='Raft L', color='blue')
line2, = ax.plot([], [], label='Raft D', color='red')

ax.set_xlim(0, len(df))
ax.set_ylim(0, df[["raft_L", "raft_D"]].values.max() * 1.1)
ax.set_title("Raft Population Growth")
ax.set_xlabel("Step")
ax.set_ylabel("Peptide Count")
ax.legend()

def update(frame):
    line1.set_data(df["step"][:frame], df["raft_L"][:frame])
    line2.set_data(df["step"][:frame], df["raft_D"][:frame])
    return line1, line2

ani = FuncAnimation(fig, update, frames=len(df), interval=20, blit=True)
ani.save(outpath, fps=30, dpi=100)
plt.close()
print(f"[✔] Raft growth animation saved to {outpath}")

