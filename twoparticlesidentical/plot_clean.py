import numpy as np
import matplotlib.pyplot as plt

# Font formatting
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

def plot_cleaned_spectrum(file_path, output_pdf_name, xlabel_str=r'$\alpha$'):
    # Load the array directly from the .npy file
    energies = np.load(file_path)  # Shape: (49, 20)
    
    # Create the alpha grid (49 points from -1.0 to 1.0)
    alphas = np.linspace(-1.0, 1.0, energies.shape[0])
    
    # Flatten/cap the spikes above 5.0
    energies_clean = np.minimum(energies, 5.0)
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
    ax.plot(alphas, energies_clean)
    
    ax.set_xlabel(xlabel_str, fontsize=14)
    ax.set_ylabel(r'Energy / $\pi^2$', fontsize=14)
    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(0.8, 5.2)

    ax.tick_params(axis='both', which='both', labelsize=14)
    ax.set_xlabel(xlabel_str, fontsize=18)
    ax.set_ylabel(r"Energy / $\pi^2$", fontsize=18)
    ax.yaxis.label.set_rotation(0)
    ax.yaxis.set_label_coords(-0.14, 0.5)
    
    plt.savefig(output_pdf_name, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_pdf_name}")


# 1. Clean plot for all alphas
plot_cleaned_spectrum(
    file_path='Energylevels_Cross_all_alphas_data.npy', 
    output_pdf_name='Energylevels_Cross_all_alphas.pdf',
    xlabel_str=r'$\alpha$'
)

# 2. Clean plot for single alpha
plot_cleaned_spectrum(
    file_path='Energylevels_Cross_single_alpha_data.npy', 
    output_pdf_name='Energylevels_Cross_single_alpha.pdf',
    xlabel_str=r'$\alpha_1$'
)