import numpy as np
import matplotlib.pyplot as plt


def load_wavefunction(filename):
    t = None
    omega = None
    channel = None

    data = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            # --- parse header
            if line.startswith("#"):
                if "t =" in line:
                    t = float(line.split("=")[1])
                elif "omega" in line:
                    omega = float(line.split("=")[1])
                elif "channel" in line:
                    channel = int(line.split("=")[1])
                continue

            # --- parse data
            parts = line.split()
            if len(parts) == 5:
                data.append([float(x) for x in parts])

    data = np.array(data)

    x = data[:, 0]
    re = data[:, 1]
    im = data[:, 2]
    density = data[:, 3]
    phase = data[:, 4]

    psi = re + 1j * im

    return x, psi, density, phase, t, omega, channel


def plot_wavefunction(x, psi, density, phase, t=None, omega=None, channel=None):
    fig, axs = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    # --- amplitude
    axs[0].plot(x, density, label=r"$|\psi|^2$")
    axs[0].plot(x, psi.real, "--", label="Re(ψ)")
    axs[0].plot(x, psi.imag, "--", label="Im(ψ)")
    axs[0].set_ylabel("Amplitude")
    axs[0].legend()
    axs[0].grid()

    # --- phase
    axs[1].plot(x, phase, label="arg(ψ)")
    axs[1].set_xlabel("x")
    axs[1].set_ylabel("Phase")
    axs[1].grid()

    title = "Wavefunction"
    if t is not None:
        title += f" | t = {t:.3f}"
    if omega is not None:
        title += f" | ω = {omega:.3f}"
    if channel is not None:
        title += f" | ch = {channel}"

    fig.suptitle(title)

    plt.tight_layout()
    plt.show()


# --- usage
if __name__ == "__main__":
    import sys

    filename = sys.argv[1] if len(sys.argv) > 1 else "wf.dat"

    x, psi, density, phase, t, omega, channel = load_wavefunction(filename)
    plot_wavefunction(x, psi, density, phase, t, omega, channel)
